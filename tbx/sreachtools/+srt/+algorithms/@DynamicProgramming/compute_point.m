function result = compute_point(obj, problem, sys, x0, varargin)

obj.validate_dependencies();
p = inputParser;
addRequired(p, 'problem', @validateproblem);
addRequired(p, 'sys', @validatesystem);
addRequired(p, 'x0');
parse(p, problem, sys, x0, varargin{:});

% 1. Ensure that the system is a Gaussian-perturbed LtiSystem
% 2. Ensure that the state dim, input_dim <=4
% 3. Input space is an axis-aligned hypercuboid
% 4. safety_tube is appropriate
% 5. optional arguments (in case of prob_str = 'first') is appropriate
otherInputHandling(sys, safety_tube, prob_str, varargin);

%% Hidden options for dynamic programming (TODO: create SReachDynOptions)
% Switch between two memory usage options
% low  - Recompute transition probability from every state to all
%        possible input and next state combinations | Low memory usage
% high - Compute all the transition probabilities at once => recursions
%        are insanely fast | High memory usage
memoryusage = 'low';
% Verbosity
verbose = 0;
% Update status every v_freq indices
v_freq = 10;

%% Compute the grid covering the target tube
% n_targets is time_horizon + 1
n_targets = length(safety_tube);

% Compute corners for gridding ==> Get corners of the largest target set in
% safety_tube
outerApproxVertices_target_sets = [];
for itt = 1:n_targets
    outerApproxVertices_target_sets = [outerApproxVertices_target_sets;
        safety_tube(itt).outerApprox.V];
end
xmax = max(outerApproxVertices_target_sets);
xmin = min(outerApproxVertices_target_sets);

% Grid computation via allcomb
if sys.state_dim == 1
    x1vec = xmin(1):x_inc:xmax(1);
    grid_x = x1vec';
    cell_xvec = {x1vec};
elseif sys.state_dim == 2
    x1vec = xmin(1):x_inc:xmax(1);
    x2vec = xmin(2):x_inc:xmax(2);
    grid_x = allcomb(x1vec,x2vec);
    cell_xvec = {x1vec,x2vec};
elseif sys.state_dim == 3
    x1vec = xmin(1):x_inc:xmax(1);
    x2vec = xmin(2):x_inc:xmax(2);
    x3vec = xmin(3):x_inc:xmax(3);
    grid_x = allcomb(x1vec,x2vec,x3vec);
    cell_xvec = {x1vec,x2vec,x3vec};
elseif sys.state_dim == 4
    x1vec = xmin(1):x_inc:xmax(1);
    x2vec = xmin(2):x_inc:xmax(2);
    x3vec = xmin(3):x_inc:xmax(3);
    x4vec = xmin(4):x_inc:xmax(4);
    grid_x = allcomb(x1vec,x2vec,x3vec,x4vec);
    cell_xvec = {x1vec,x2vec,x3vec,x4vec};
end
% No. of grid points
n_grid_x = length(grid_x);

%% Input gridding
umax = max(sys.input_space.V);
umin = min(sys.input_space.V);
if sys.input_dim == 1
    grid_u = allcomb(umin(1):u_inc:umax(1));
elseif sys.input_dim == 2
    grid_u = allcomb(umin(1):u_inc:umax(1), ...
                     umin(2):u_inc:umax(2));
elseif sys.input_dim == 3
    grid_u = allcomb(umin(1):u_inc:umax(1), ...
                     umin(2):u_inc:umax(2), ...
                     umin(3):u_inc:umax(3));
elseif sys.input_dim == 4
    grid_u = allcomb(umin(1):u_inc:umax(1), ...
                     umin(2):u_inc:umax(2), ...
                     umin(3):u_inc:umax(3), ...
                     umin(4):u_inc:umax(4));
end

%% For trapezoid rule, we penalize 1/2 per dimension at the endpoints
% Endpoint iff one of the dimensions is xmin or xmax
% Use max to implement OR and sum to count how many active dimensions
n_active_dims = max(sum(grid_x == xmin,2),sum(grid_x == xmax,2));
% Compute scaling for the trapezoidal rule
fraction_at_grid = 2.^(-n_active_dims);
% Create dx with appropriate scaling based on where it is located
delta_x_grid = (x_inc^sys.state_dim).* fraction_at_grid;


%% Initialize a matrix to store the value functions
mat_prob_x = zeros(n_targets, n_grid_x);

%% How do we compute the transition probability
switch memoryusage
    case 'high'
        % Compute transition probabilities for every x, u combination
        transition_prob_with_delta_all = computeTransProbWithDelta(sys, ...
            grid_x, grid_u, delta_x_grid, verbose, v_freq);
        % Return a n_grid_x * n_grid_u matrix corresponding current state
        % grid_x(ix)
        transition_prob_with_delta = @(ix) ...
            transition_prob_with_delta_all{ix}';
    case 'low'
        % Compute the transition probability at every step
        transition_prob_with_delta = @(ix) ...
            computeTransProbWithDeltaAtX(sys, ix, grid_x, grid_u, ...
                delta_x_grid)';
end

%% Implement dynamic programming
switch prob_str
    case 'term'
        terminal_indicator_x = safety_tube(n_targets).contains(grid_x');
        % fprintf('Set optimal value function at t=%d\n',n_targets-1);
        mat_prob_x(n_targets,:) = terminal_indicator_x;

        for itt = n_targets - 1:-1:1
            % fprintf('Compute optimal value function at t=%d\n', itt - 1);
            % Obtain V_{t+1}
            old_prob_x = mat_prob_x(itt+1,:);
            % Check which of the grid points need to be iterated over
            current_indicator_x = safety_tube(itt).contains(grid_x');
            % Verbosity: Initial display
            if verbose
                switch memoryusage
                    case 'low'
                        n_ix = nnz(current_indicator_x);
                        fprintf(['Time t = %d | Computing V_t(x)....', ...
                            '%3.2f%%'], itt - 1, round(0/n_ix*100,2));
                    case 'high'
                        fprintf('Recursion at t = %d\n', itt - 1);
                end
            end
            % Iterate over all these points and compute
            % max_u \int_X V_{t+1}(x_{t+1}(u))
            %                       * transitionProb(x_{t+1},u)dx_{t+1}
            for ix = find(current_indicator_x==1)
                mat_prob_x(itt,ix) = max(...
                    old_prob_x*transition_prob_with_delta(ix));
                % Verbosity: continue_display
                if verbose && strcmpi(memoryusage,'low') && ~mod(ix, v_freq)
                    percent_completed = round(ix/n_ix*100,2);
                    if percent_completed < 10
                        fprintf('\b\b\b\b\b%3.2f%%', percent_completed);
                    else
                        fprintf('\b\b\b\b\b\b%3.2f%%', percent_completed);
                    end
                end
            end
            if verbose && strcmpi(memoryusage,'low')
                fprintf('\b\b\b\b\b\b%3.2f%%\n', 100);
            end
        end
    case 'first'
        target_set = varargin{1};
        target_set_indicator_x = target_set.contains(grid_x');
        indx_grid_points_outside_target_set = ~double(target_set_indicator_x);

        % Initialize
        % fprintf('Set optimal value function at t=%d\n',n_targets-1);
        mat_prob_x(n_targets,:) = target_set_indicator_x;

        for itt = n_targets - 1:-1:1
            % fprintf('Compute optimal value function at t=%d\n', itt - 1);
            % Obtain V_{t+1}
            old_prob_x = mat_prob_x(itt+1,:);
            % Set points that reached the target set as one
            mat_prob_x(itt,:) = target_set_indicator_x;
            % Iterate over the remaining safe points and compute
            % max_u \int_X V_{t+1}(x_{t+1}(u))transitionProb(x_{t+1},u)dx_{t+1}
            current_indicator_x = indx_grid_points_outside_target_set.* ...
                safety_tube(itt).contains(grid_x');
            for ix = find(current_indicator_x==1)
                mat_prob_x(itt,ix) = max(old_prob_x * ...
                    transition_prob_with_delta{ix}');
            end
        end
end
prob_x = mat_prob_x(1,:);
varargout{1} = cell_xvec;
varargout{2} = grid_x;
varargout{3} = mat_prob_x;

end
