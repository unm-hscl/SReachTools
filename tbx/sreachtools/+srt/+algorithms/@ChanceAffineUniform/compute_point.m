function results = compute_point(obj, prob, sys, x0, varargin)
% COMPUTE_POINT Computes a point solution for the ChanceAffineUniform algorithm.
%
%   results = COMPUTE_POINT(obj, problem, sys, x0)
%   results = COMPUTE_POINT(obj, problem, sys, x0, ...)
%
%   [lb_stoch_reach, opt_input_vec, opt_input_gain, risk_alloc_state, ...
%       risk_alloc_input] = SReachPointCcAu(sys, initial_state, safety_tube, ...
%       options)
%
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.StateDimension)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - Collection of user-specified options for 'chance-affine-uni'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   lb_stoch_reach
%               - Lower bound on the stochastic reachability of a target tube
%                 problem computed using chance constraints and
%                 difference-of-convex techniques
%   opt_input_vec,
%     opt_input_gain
%               - Controller U=MW+d for a concatenated input vector
%                   U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
%                   vector W=[w_0; w_1; ...; w_{N-1}].
%                   - opt_input_gain: Affine controller gain matrix of dimension
%                       (sys.InputDimension*N) x (sys.DisturbanceDimension*N)
%                   - opt_input_vec: Open-loop controller: column vector
%                     dimension
%                       (sys.InputDimension*N) x 1
%   risk_alloc_state
%               - Risk allocation for the state constraints
%   risk_alloc_input
%               - Risk allocation for the input constraints
%
% See also SReachPoint.
%
% Notes:
% * See @LtiSystem/getConcatMats for more information about the notation used.
%
% ============================================================================
%
%
%
%   This algorithm is part of the Stochastic Reachability Toolbox.
%   License for the use of this algorithm is given in
%   https://sreachtools.github.io/license/

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

% Get half space representation of the target tube and time horizon
% skipping the first time step
[concat_safety_tube_A, concat_safety_tube_b] = safety_tube.concatenate(...
    [2 N+1]);

%% Halfspace-representation of U^N, H, G,mean_X_sans_input, cov_X_sans_input
% GUARANTEES: Non-empty input sets (polyhedron)
[concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, ...
    N);

% GUARANTEES: Compute the input concat and disturb concat transformations
[~, H, G] = getConcatMats(sys, N);

% GUARANTEES: well-defined initial_state and N
sysnoi = srt.systems.LtvSystem( ...
    'A', @(t) sys.A(t), ...
    'F', sys.F, ...
    'w', sys.Disturbance);

X_sans_input_rv = srt.SReachFwd('concat-stoch', sysnoi, ...
    x0, N);
mean_X_zi = X_sans_input_rv.Mean();
mean_X_zi = mean_X_zi(sysnoi.StateDimension + 1:end);

mean_W = kron(ones(N,1), sys.Disturbance.Mean());

%% Compute M --- the number of polytopic halfspaces to worry about
n_lin_state = size(concat_safety_tube_A,1);
n_lin_input = size(concat_input_space_A,1);

%% Covariance of W vector
cov_concat_disturb = kron(eye(N), sys.Disturbance.Sigma());

% Compute a sparse square root of a matrix: chol produces a sqrt such that
% sqrt' * sqrt = M. Hence, whenever, we left multiply (inside a norm), we
% must transpose.
[sqrt_cov_concat_disturb, p] = chol(cov_concat_disturb);
if p > 0
    % Non-positive definite matrix can not use Cholesky's decomposition
    % Use sqrt to obtain a symmetric non-sparse square-root matrix
    sqrt_cov_concat_disturb = sqrt(cov_concat_disturb);
end

lb_stoch_reach = -1;
opt_input_vec = nan(sys.InputDimension * N,1);
opt_input_gain = [];
risk_alloc_state = nan(n_lin_state,1);
risk_alloc_input = nan(n_lin_input,1);

%% Bisect our ways into the risk allocation
state_prob_lb = 0;
state_prob_ub = 1;
while state_prob_ub - state_prob_lb > obj.state_bisect_tol
    % Compute the risk allocations for all of the state constraints
    state_prob_test = (state_prob_ub + state_prob_lb)/2;
    state_viol_risk_per_ineq = (1 - state_prob_test)/n_lin_state;

    % Within a specified state constraint risk of violation, we search for
    % a controller that performs satisfactorily up to a maximum allowed
    % input constraint risk of violation. If possible, we seek for a
    % smaller input_viol_prob_lb since that increases our lower bound
    % on the reach probability
    input_viol_prob_lb = 0;
    % Ensures that Theorem 1 hypothesis is satisfied
    input_viol_prob_ub = min(obj.max_input_viol_prob, state_prob_test);

    % Keep track, if we obtained a feasible solution among all possible
    % input constraint risk violators => we can dream higher in terms
    % of state constraint probability
    at_least_once_feasible = 0;
    while input_viol_prob_ub - input_viol_prob_lb > obj.input_bisect_tol
        % Compute the risk allocations for all of the input constraints
        input_viol_prob_test = (input_viol_prob_lb + input_viol_prob_ub)/2;
        input_viol_risk_per_ineq = input_viol_prob_test/n_lin_input;

        % The iteration values are updated at the end of the problem
        cvx_begin quiet
            variable M(sys.InputDimension*N,sys.DisturbanceDimension*N)
            variable d(sys.InputDimension * N, 1);
            maximize 0;
            subject to
                % Causality constraints on M_matrix
                for time_indx = 1:N
                    M((time_indx-1)*sys.InputDimension + 1:...
                        time_indx*sys.InputDimension, ...
                        (time_indx-1)*sys.DisturbanceDimension+1:end) == 0;
                end
                % Individual chance constraint --- input
                concat_input_space_A *  (M * mean_W + d) + ...
                    norminv(1-input_viol_risk_per_ineq) ...
                    * norms(concat_input_space_A* M * ...
                        sqrt_cov_concat_disturb',2,2)<=concat_input_space_b;
                % Individual chance constraint --- state
                % Mean trajectory constraint (mean_X_zi already has Zx + GW)
                concat_safety_tube_A * (mean_X_zi + H * (M * mean_W + d))...
                    + norminv(1-state_viol_risk_per_ineq)...
                    * norms(concat_safety_tube_A* (H * M + G) * ...
                        sqrt_cov_concat_disturb',2,2)<= concat_safety_tube_b
        cvx_end

        obj.print_verbose(1, ...
            ['Safety prob: %1.3f, Input viol: %1.3f | ', ...
            'CVX status: %s\n'], state_prob_test, ...
            input_viol_prob_test, cvx_status);

        switch cvx_status
            case 'Solved'
                % Raise the expected safety probability
                state_prob_lb = state_prob_test;
                results.lb_stoch_reach = (state_prob_lb - input_viol_prob_test)./...
                    (1 - input_viol_prob_test);
                results.risk_alloc_state = state_viol_risk_per_ineq * ...
                    ones(n_lin_state,1);
                results.risk_alloc_input = input_viol_risk_per_ineq * ...
                    ones(n_lin_input,1);
                results.opt_input_vec = d;
                results.opt_input_gain = M;
                at_least_once_feasible = 1;
                break;
            otherwise
                % Allow for more input_violations
                input_viol_prob_lb = input_viol_prob_test;
        end
    end
    if ~at_least_once_feasible
        % Lower the expected safety probability
        state_prob_ub = state_prob_test;
    end
end

end
