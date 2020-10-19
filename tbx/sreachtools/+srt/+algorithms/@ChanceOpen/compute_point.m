function results = compute_point(obj, prob, sys, x0, varargin)
% Solve the problem of stochastic reachability of a target tube (a lower bound
% on the maximal reach probability and an open-loop controller synthesis) using
% convex chance-constrained optimization
% =============================================================================
%
% SReachPointCcO implements a chance-constrained convex underapproximation to
% the stochastic reachability of a target tube prob. The original problem was
% formulated (for the simpler problem of terminal hitting-time stochastic
% reach-avoid problem) in
%
% K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
% spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
% 2013.
%
% This function implements a convex solver-friendly using piecewise-affine
% overapproximations of the convex constraints, as discussed in
%
% A. Vinod and M. Oishi. Affine controller synthesis for stochastic reachability
% via difference of convex programming. In Proc. Conf. Dec. & Ctrl., 2019.
% (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
%
%    High-level desc.   : Use Boole's inequality, Gaussian random vector, and
%                         piecewise linear approximation of the inverse of the
%                         standard normal cumulative density function to create
%                         a linear program-based approximation to the original
%                         optimization
%    Approximation      : Guaranteed underapproximation
%    Controller type    : Open-loop controller that satisfies the hard
%                         input bounds
%    Optimality         : Optimal open-loop controller for the
%                         underapproximation problem due to convexity guarantees
%
% =============================================================================
%
% [lb_stoch_reach, opt_input_vec, risk_alloc_state, varargout] = ...
%    SReachPointCcO(sys, initial_state, safety_tube, options)
%
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - Collection of user-specified options for 'chance-open'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   lb_stoch_reach
%               - Lower bound on the stochastic reachability of a target tube
%                 problem computed using convex chance constraints and
%                 piecewise affine approximation. While it is expected to lie in
%                 [0,1], it is set to -1 in cases where the CVX optimization
%                 fails (cvx_status \neq Solved).
%   opt_input_vec
%               - Open-loop controller: column vector of dimension
%                 (sys.input_dim*N) x 1
%   risk_alloc_state
%               - Risk allocation for the state constraints
%   extra_info  - [Optional] Useful information to construct the
%                   reachability problem | Used by 'genzps-open' to avoid
%                   unnecessary recomputation
%                 Matlab struct with members --- concat_safety_tube_A,
%                   concat_safety_tube_b, concat_input_space_A,
%                   concat_input_space_b, H, mean_X_sans_input,
%                   cov_X_sans_input;
%
% See also SReachPoint.
%
% Notes:
% * We recommend using this function through SReachPoint.
% * This function requires CVX to work.
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

% Target tubes has polyhedra T_0, T_1, ..., T_{N}
N = prob.TimeHorizon - 1;


% Get half space representation of the target tube and time horizon
% skipping the first time step
[concat_safety_tube_A, concat_safety_tube_b] = ...
    prob.TargetTube.concatenate([2 N+1]);

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

X_sans_input_rv = srt.SReachFwd('concat-stoch', sysnoi, x0, N);

mean_X_sans_input = X_sans_input_rv.Mean();
mean_X_sans_input = mean_X_sans_input(sysnoi.StateDimension+1:end);

cov_X_sans_input = X_sans_input_rv.Sigma();
cov_X_sans_input = cov_X_sans_input( ...
    sysnoi.StateDimension+1:end, ...
    sysnoi.StateDimension+1:end);


%% Compute M --- the number of polytopic halfspaces to worry about
n_lin_state = size(concat_safety_tube_A,1);

%% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i} = ||\sqrt\Sigma_X*h_i||
% cholesky > cov_X_sans_input = sqrt_cov_X_sans_input'*sqrt_cov_X_sans_input
[sqrt_cov_X_sans_input, p] = chol(cov_X_sans_input);
if p > 0
    % Non-positive definite matrix can not use Cholesky's decomposition
    % Use sqrt to obtain a symmetric non-sparse square-root matrix
    sqrt_cov_X_sans_input = sqrt(cov_X_sans_input);
end
% Hence, we need the transpose on sqrt_cov_X
scaled_sigma_vec = norms(concat_safety_tube_A*sqrt_cov_X_sans_input', 2, 2);


%% Obtain the piecewise linear overapproximation of norminvcdf in [0,0.5]
[invcdf_approx_m, invcdf_approx_c, lb_deltai] =...
    computeNormCdfInvOverApprox(0.5, obj.pwa_accuracy, n_lin_state);

%% Solve the feasibility problem
cvx_begin quiet
    variable U_vector(sys.InputDimension * N, 1);
    variable mean_X(sys.StateDimension * N, 1);
    variable deltai(n_lin_state, 1);
    variable norminvover(n_lin_state, 1);
    minimize sum(deltai)
    subject to
        mean_X == mean_X_sans_input + H * U_vector;
        concat_input_space_A * U_vector <= concat_input_space_b;
        for deltai_indx = 1:n_lin_state
            norminvover(deltai_indx) >= invcdf_approx_m.* ...
                deltai(deltai_indx) + invcdf_approx_c;
        end
        concat_safety_tube_A * mean_X + scaled_sigma_vec.* norminvover ...
            <= concat_safety_tube_b;
        deltai >= lb_deltai;
        deltai <= 0.5;
cvx_end

if strcmpi(cvx_status, 'Solved')
    results.lb_stoch_reach = 1-sum(deltai);
    results.opt_input_vec = U_vector;
    results.risk_alloc_state = deltai;
else
    results.lb_stoch_reach = -1;
    results.opt_input_vec = nan(sys.InputDimension * N, 1);
    results.risk_alloc_state = nan(n_lin_state, 1);
end

%% Create the other info for use if necessary
results.extra_info.concat_safety_tube_A = concat_safety_tube_A;
results.extra_info.concat_safety_tube_b = concat_safety_tube_b;
results.extra_info.concat_input_space_A = concat_input_space_A;
results.extra_info.concat_input_space_b = concat_input_space_b;
results.extra_info.H = H;
results.extra_info.mean_X_sans_input = mean_X_sans_input;
results.extra_info.cov_X_sans_input = cov_X_sans_input;

end
