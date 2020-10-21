function results = compute_point(obj, prob, sys, x0, varargin)
% Solve the problem of stochastic reachability of a target tube (a lower bound
% on the maximal reach probability and an open-loop controller synthesis) using
% particle filter control
% =============================================================================
%
% SReachPointPaO implements a mixed-integer linear program-based approximation
% to the stochastic reachability of a target tube prob. This solution is
% based off the particle filter control formulation (for the simpler terminal
% hitting-time stochastic reach-avoid problem) discussed in
%
% K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
% spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
% 2013.
%
%    High-level desc.   : Sample scenarios based on the additive noise and solve
%                         a mixed-integer linear program to make the maximum
%                         number of scenarios satisfy the reachability objective
%    Approximation      : No direct approximation guarantees. Accuracy improves
%                         as the number of scenarios considered increases.
%    Controller type    : Open-loop controller that satisfies the hard input
%                         bounds
%    Optimality         : Optimal (w.r.t scenarios drawn) open-loop controller
%                         for the underapproximation problem
%
% =============================================================================
%
% [approx_stoch_reach, opt_input_vec] = SReachPointPaO(sys, x0, ...
%   safety_tube, options)
%
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   x0- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.StateDimension)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - Collection of user-specified options for 'particle-open'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   approx_stoch_reach
%               - An approximation of the stochastic reachability of a target
%                 tube problem computed using particle control. While it is
%                 expected to lie in [0,1], it is set to -1 in cases where the
%                 CVX optimization fails (cvx_status \not\in {Solved,
%                 Inaccurate/Solved}).
%   opt_input_vec
%               - Open-loop controller: column vector of dimension
%                 (sys.InputDimension*N) x 1
%
% See also SReachPoint.
%
% Notes:
% * This function requires CVX with Gurobi as the backend solver for optimizing
%   the resulting mixed-integer linear program.
% * See @LtiSystem/getConcatMats for more information about the notation used.
%
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%

p = inputParser;
addRequired(p, 'prob', @obj.validateproblem);
addRequired(p, 'sys', @obj.validatesystem);
addRequired(p, 'x0');
parse(p, prob, sys, x0, varargin{:});

results = struct;

import srt.*

safety_tube = prob.TargetTube;

% Target tubes has polyhedra T_0, T_1, ..., T_{N}
N = prob.TimeHorizon - 1;

approx_stoch_reach = -1;
opt_input_vec = nan(sys.InputDimension * N,1);

% Requires Gurobi since we are solving a MILP
[default_solver, solvers_cvx] = cvx_solver;
n_Gurobi_solver = nnz(contains(solvers_cvx,'Gurobi'));

n_particles = obj.n_particles;

if n_Gurobi_solver == 0

    warning('SReachTools:runtime',['SReachPointVoO returns a trivial ', ...
        'result since Gurobi (a MILP solver) was not setup.']);

else

    if ~contains(default_solver, 'Gurobi') && n_Gurobi_solver >= 2
        warning('SReachTools:runtime', sprintf(['SReachPointVoO ', ...
            'requires a MILP solver. Found %d Gurobi solvers.\n', ...
            'Choosing the Gurobi solver bundled with CVX.\nSet the ', ...
            'desired Gurobi solver, before calling this function.'], ...
            n_Gurobi_solver));
    end


    % Get half space representation of the target tube and time horizon
    % skipping the first time step
    [concat_safety_tube_A, concat_safety_tube_b] = ...
        safety_tube.concatenate([2 N+1]);

    % Normalize the hyperplanes so that Ax <= b + M(1-bin_x) can be
    % implemented with M = 100 (numerically stable compared to using large
    % M). Here, b<=1e-3 or 1
    [concat_safety_tube_A, concat_safety_tube_b] = ...
        normalizeForParticleControl(concat_safety_tube_A, ...
            concat_safety_tube_b);

    % Halfspace-representation of U^N, and matrices Z, H, and G
    % GUARANTEES: Non-empty input sets (polyhedron)
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(...
        sys, N);
    % Compute the input concatenated transformations
    [Z, H, G] = getConcatMats(sys, N);

    % Compute M --- the number of polytopic halfspaces to worry about
    n_lin_state = size(concat_safety_tube_A,1);

    obj.print_verbose(1, 'Required number of particles: %d\n', ...
        n_particles);

    obj.print_verbose(1, 'Creating random variable realizations....');

    % Compute the stochasticity of the concatenated disturbance random vec
    W = concat(sys.Disturbance, N);
    % Create realizations of W arranged columnwise
    W_realizations = W.sample(n_particles);

    obj.print_verbose(1, 'Done\n');

    % Implementation of Problem 2 in Lesser CDC 2013
    % Solve the mixed-integer linear program
    obj.print_verbose(2, 'Objective value needs to be scaled by %1.3f\n', ...
        1/n_particles);

    obj.print_verbose(1, 'Setting up CVX prob....');

    cvx_begin quiet

        if ~contains(default_solver,'Gurobi')
            cvx_solver Gurobi;
        end
        variable U_vector(sys.InputDimension * N,1);
        variable mean_X(sys.StateDimension * N,n_particles);
        variable bin_x(1,n_particles) binary;
        maximize sum(bin_x)/n_particles
        subject to
            mean_X == repmat(Z * x0 + H * U_vector, ...
                1, n_particles) + G * W_realizations;
            concat_input_space_A * U_vector <= concat_input_space_b;
            concat_safety_tube_A * mean_X <= repmat( ...
                concat_safety_tube_b, 1, n_particles) + ...
                obj.bigM * repmat(1-bin_x,n_lin_state,1);


        obj.print_verbose(1, 'Done\nParsing and solving the MILP....');

    cvx_end

    obj.print_verbose(1, 'Done\n');

    %% Overwrite the solutions
    switch cvx_status
        case {'Solved', 'Inaccurate/Solved'}
            results.approx_stoch_reach = sum(bin_x)/n_particles;
            results.opt_input_vec = U_vector;
        otherwise

    end

end

end
