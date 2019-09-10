function result = compute_point(obj, problem, sys, x0, varargin)

p = inputParser;
addRequired(p, prb, @obj.validateproblem);
addRequired(p, sys, @obj.validatesystem);
parse(p, prb, sys);

% Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
time_horizon = length(safety_tube) - 1;

otherInputHandling(sys, options, time_horizon);

% Get half space representation of the target tube and time horizon
% skipping the first time step
[concat_safety_tube_A, concat_safety_tube_b] = safety_tube.concat(...
    [2 time_horizon+1]);

%% Halfspace-representation of U^N, H, G,mean_X_sans_input, cov_X_sans_input
% GUARANTEES: Non-empty input sets (polyhedron)
[concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, ...
    time_horizon);
% GUARANTEES: Compute the input concat and disturb concat transformations
[~, H, G] = getConcatMats(sys, time_horizon);
% GUARANTEES: well-defined initial_state and time_horizon
sysnoi = LtvSystem('StateMatrix',sys.state_mat,'DisturbanceMatrix', ...
    sys.dist_mat,'Disturbance',sys.dist);
X_sans_input_rv = SReachFwd('concat-stoch', sysnoi, initial_state, ...
    time_horizon);
mean_X_zi = X_sans_input_rv.mean();
mean_X_zi = mean_X_zi(sysnoi.state_dim + 1:end);

mean_W = kron(ones(time_horizon,1), sys.dist.mean());


%% Compute M --- the number of polytopic halfspaces to worry about
n_lin_state = size(concat_safety_tube_A,1);
n_lin_input = size(concat_input_space_A,1);

%% Covariance of W vector
cov_concat_disturb = kron(eye(time_horizon),sys.dist.cov());
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
opt_input_vec = nan(sys.input_dim * time_horizon,1);
opt_input_gain = [];
risk_alloc_state = nan(n_lin_state,1);
risk_alloc_input = nan(n_lin_input,1);

%% Bisect our ways into the risk allocation
state_prob_lb = 0;
state_prob_ub = 1;
while state_prob_ub - state_prob_lb > options.state_bisect_tol
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
    input_viol_prob_ub = min(options.max_input_viol_prob, state_prob_test);

    % Keep track, if we obtained a feasible solution among all possible
    % input constraint risk violators => we can dream higher in terms
    % of state constraint probability
    at_least_once_feasible = 0;
    while input_viol_prob_ub - input_viol_prob_lb > options.input_bisect_tol
        % Compute the risk allocations for all of the input constraints
        input_viol_prob_test = (input_viol_prob_lb + input_viol_prob_ub)/2;
        input_viol_risk_per_ineq = input_viol_prob_test/n_lin_input;

        % The iteration values are updated at the end of the problem
        cvx_begin quiet
            variable M(sys.input_dim*time_horizon,sys.dist_dim*time_horizon)
            variable d(sys.input_dim * time_horizon, 1);
            maximize 0;
            subject to
                % Causality constraints on M_matrix
                for time_indx = 1:time_horizon
                    M((time_indx-1)*sys.input_dim + 1:...
                        time_indx*sys.input_dim, ...
                        (time_indx-1)*sys.dist_dim+1:end) == 0;
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
        if options.verbose >= 1
            fprintf(['Safety prob: %1.3f, Input viol: %1.3f | ', ...
                'CVX status: %s\n'], state_prob_test, ...
                input_viol_prob_test, cvx_status);
        end
        switch cvx_status
            case 'Solved'
                % Raise the expected safety probability
                state_prob_lb = state_prob_test;
                lb_stoch_reach = (state_prob_lb - input_viol_prob_test)./...
                    (1 - input_viol_prob_test);
                risk_alloc_state = state_viol_risk_per_ineq * ...
                    ones(n_lin_state,1);
                risk_alloc_input = input_viol_risk_per_ineq * ...
                    ones(n_lin_input,1);
                opt_input_vec = d;
                opt_input_gain = M;
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
