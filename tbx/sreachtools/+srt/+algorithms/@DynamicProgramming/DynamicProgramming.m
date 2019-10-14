classdef DynamicProgramming < srt.algorithms.Algorithm
% DYNAMICPROGRAMMING Dynamic programming implementation.

    properties
        % State space grid.
        ss_grid_ double {mustBeNumeric}
        
        % Input space grid.
        is_grid_
        
        % Probability values at state space grid points.
        prob_ double {mustBeNumeric}

        ss_grid_is_rect_ = []

        is_grid_is_rect_ = []

        int_dx_ = []
    end

    properties (Access = private, Dependent)
        xmax_

        xmin_
    end

    properties (Dependent)
        % STATESPACEGRID State space grid.
        StateSpaceGrid
        
        % INPUTSPACEGRID Input space grid.
        InputSpaceGrid

        Probability
    end

    methods
        function obj = DynamicProgramming(varargin)
            % DYNAMICPROGRAMMING Construct an instance of the algorithm.

            p = inputParser();
            addParameter(p, 'StateSpaceGrid', []);
            addParameter(p, 'InputSpaceGrid', []);
            addParameter(p, 'verbose', 0, @(x) validateattributes(x, ...
                {'numeric'}, {'scalar', 'integer', 'nonnegative'}));
            addParameter(p, 'MemoryUsage', 'High', @(x) validatestring(x, ...
                {'Low', 'High'}));
            addParameter(p, 'NonGaussianParticleCount', 10000, ...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'nonnegative'}));

            parse(p, varargin{:});

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            obj.StateSpaceGrid = p.Results.StateSpaceGrid;
            obj.InputSpaceGrid = p.Results.InputSpaceGrid;

            
        end
    end

    methods
        function val = get.xmax_(obj)
            val = max(obj.StateSpaceGrid);
        end

        function val = get.xmin_(obj)
            val = min(obj.StateSpaceGrid);
        end

        function grid = get.StateSpaceGrid(obj)
            grid = obj.ss_grid_;
        end

        function set.StateSpaceGrid(obj, grid)
            obj.ss_grid_ = obj.check_and_set_grid_(grid);
            obj.ss_grid_is_rect_ = obj.grid_is_rect_(obj.ss_grid_);
        end

        function grid = get.InputSpaceGrid(obj)
            grid = obj.is_grid_;
        end
        

        function set.InputSpaceGrid(obj, grid)
            obj.is_grid_ = obj.check_and_set_grid_(grid);
            obj.is_grid_is_rect_ = obj.grid_is_rect_(obj.is_grid_);
        end

        function result = compute_point(obj, prob, sys, ~)
            obj.validatesystem(sys);
            obj.validateproblem(prob);

            N = prob.TimeHorizon;
            result = struct();
            if isempty(obj.is_grid_);
                result.policy = [];
            else
                result.policy = zeros(size(obj.ss_grid_, 1), N-1);
            end

            if isa(prob, 'srt.problems.FirstHitting')
                % First hitting time recursion
                result.prob = zeros(size(obj.ss_grid_, 1), N);
                obj.print_verbose(1, 'Computing V_%d(x).\n', N);
                result.prob(:, N) = contains(prob.TargetTube(N), obj.ss_grid_');

                for lv = N-1:-1:1
                    % Recursion
                    obj.print_verbose(1, 'Computing V_%d(x).\n', lv);
                    safe_set = prob.ConstraintTube(lv);
                    target_set = prob.TargetTube(lv);
                    for ix = 1:size(obj.ss_grid_, 1)
                        if target_set.contains(obj.ss_grid_(ix, :)')
                            continue;
                        end
                        
                        if ~safe_set.contains(obj.ss_grid_(ix, :)')
                            continue;
                        end

                        if isempty(obj.is_grid_)
                            xp = sys.onestep(lv, obj.ss_grid_(ix, :)');
                            tx = obj.grid_transition_prob(xp);
                            result.prob(ix, lv) = tx' * result.prob(:, lv+1);
                        else
                            for iu = 1:size(obj.is_grid_, 1)
                                xp = sys.onestep(lv, obj.ss_grid_(ix, :)', ...
                                    obj.is_grid_(iu, :)');
                                
                                tx = obj.grid_transition_prob(xp);
                                p = tx' * result.prob(:, lv+1);
                                if p > result.prob(ix, lv)
                                    result.policy(ix, lv) = iu;
                                    result.prob(ix, lv) = p;
                                end
                            end
                        end
                    end
            elseif isa(prob, 'srt.problems.TerminalHitting')
                % Terminal hitting time recursion
                result.prob = zeros(size(obj.ss_grid_, 1), N);
                obj.print_verbose(1, 'Computing V_%d(x).\n', N);
                result.prob(:, N) = contains(prob.TargetTube(N), obj.ss_grid_');

                for lv = N-1:-1:1
                    % Recursion
                    obj.print_verbose(1, 'Computing V_%d(x).\n', lv);
                    safe_set = prob.ConstraintTube(lv);
                    for ix = 1:size(obj.ss_grid_, 1)
                        if ~safe_set.contains(obj.ss_grid_(ix, :)')
                            continue;
                        end

                        if isempty(obj.is_grid_)
                            xp = sys.onestep(lv, obj.ss_grid_(ix, :)');
                            tx = obj.grid_transition_prob(xp);
                            result.prob(ix, lv) = tx' * result.prob(:, lv+1);
                        else
                            for iu = 1:size(obj.is_grid_, 1)
                                xp = sys.onestep(lv, obj.ss_grid_(ix, :)', ...
                                    obj.is_grid_(iu, :)');
                                
                                tx = obj.grid_transition_prob(xp);
                                p = tx' * result.prob(:, lv+1);
                                if p > result.prob(ix, lv)
                                    result.policy(ix, lv) = iu;
                                    result.prob(ix, lv) = p;
                                end
                            end
                        end
                    end
                end

            end
        end

        function tx = grid_transition_prob(obj, rv)
            persistent nongauss_warning_raised
            if isempty(nongauss_warning_raised)
                nongauss_warning_raised = false;
            end

            if obj.ss_grid_is_rect_ && ...
                isa(rv, 'srt.disturbances.Gaussian')
                % Gaussian disturbances can be computed faster by an 
                % implementation of trapezoidal integration
                
                if isempty(obj.int_dx_)
                    obj.set_int_dx_();
                end
                
                tx = rv.pdf(obj.ss_grid_) .* obj.int_dx_;
                
            else
                
                % Transition probabilites for generic disturbances must be
                % computed empirically

                if ~nongauss_warning_raised
                    warning(['Dynamic programming for non-Gaussian ', ...
                        'disturbances is still under construction. Results ', ...
                        'be innacurate.']);
                    nongauss_warning_raised = true;
                end
                
                tx = zeros(size(obj.ss_grid_, 1), 1);
                samples = rv.sample(obj.NonGaussianParticleCount);

                for s = samples
                    [~, i] = min(sum(abs(obj.ss_grid_' - s)));
                    tx(i) = tx(i) + 1;
                end

                tx = tx / obj.NonGaussianParticleCount;

            end
        end

        % function transition_prob_with_delta_at_x = computeTransProbWithDeltaAtX(sys, ...
        %     ix, grid_x, grid_u, delta_x_grid)
        
        %     % Get transition probability for all points
        %     transition_prob_with_delta_at_x = zeros(length(grid_u), length(grid_x));
            
        %     % Compute the random vector (x_{k+1} - B u_k), given by A x_k + F w_k
        %     next_x_zi = sys.state_mat * grid_x(ix,:)' + sys.dist_mat * sys.dist;   
            
        %     for iu = 1:length(grid_u)
        %         % Compute the random vector x_{k+1} under a specific u_k
        %         next_x = next_x_zi + sys.input_mat * grid_u(iu,:);
        %         % Transition probability is the pdf times delta for integration
        %         transition_prob_with_delta_at_x(iu,:) = mvnpdf(grid_x, ...
        %             next_x.mean()', next_x.cov())' .* delta_x_grid';
        %     end
        % end
    end

    methods (Access = private)
        function validateproblem(obj, prob)
            validateattributes(prob, {'srt.problems.FirstHitting', ...
                'srt.problems.TerminalHitting'}, {'nonempty'});
        end

        function g = check_and_set_grid_(obj, grid)
            if isempty(grid)
                g = grid;
                return;
            end

            if isnumeric(grid) && ismatrix(grid)
                if any(size(grid) == 1)
                    g = reshape(grid, [], 1);
                else
                    g = grid;
                end

                return;
            end

            if iscell(grid)
                % Three ways of specifying the grid using a cell array
                %   1) {integer=n, vector}
                %       Builds with [X1, ..., Xn] = ndgrid(vector)
                %   2) {vector1, vector2, ...}
                %       Builds with [X1, X2, ...] = ndgrid(vector1, vector2, ...)
                %   3) {X1, X2, ...}
                %       Assumes X1, X2, ... are formed from meshgrid or ndgrid.
                %       Just reshapes
                % 

                valid_grid_input = true;

                if (isnumeric(grid{1}) && isscalar(grid{1})) && ...
                    (length(grid) == 2 && isnumeric(grid{2}) ...
                        && isvector(grid{2}))

                    scalar_vector = true;
                else
                    scalar_vector = false;
                end

                if ~scalar_vector
                    all_vectors = true;
                    all_grids   = true;

                    for sub = grid
                        if isnumeric(sub) && isvector(sub)
                            all_vectors = all_vectors * true;
                            all_grids   = all_grids * false;
                        elseif isnumeric(sub) && length(size(sub)) == length(grid)
                            all_vectors = all_vectors * false;
                            all_grids   = all_grids * true;
                        else
                            all_vectors = all_vectors * false;
                            all_grids   = all_grids * false;
                            break;
                        end
                    end
                end

                if scalar_vector
                    M = cell(1, grid{1});
                    [M{:}] = ndgrid(grid{2});
                elseif all_vectors
                    M = cell(1, grid{1});
                    [M{:}] = ndgrid(grid{:});
                elseif all_grids
                    M = grid;
                else
                    error('Invalid grid input: Please see help.');
                end

                g = zeros(numel(M{1}), length(M));
                for i = 1:length(M)
                    g(:, i) = reshape(M{i}, [], 1);
                end
            end
        end

        function bl = grid_is_rect_(obj, g)
            if isempty(g)
                bl = true;
            end

            bl = false;
            if floor(sqrt(size(g, 1))) == ceil(sqrt(size(g, 1)))
                rect_dim_count = 0;
                for lv = 1:size(g, 2)
                    rg = reshape(g(:, lv), sqrt(size(g(:, lv), 1)) * ones(1, 2));
                    drg1 = diff(rg, 1, 1);
                    drg2 = diff(rg, 1, 2);

                    if all(abs(drg1 - drg1(1) < 1e-6), 'all') && ...
                        all(abs(drg2 - drg2(1)) < 1e-6, 'all')
                        rect_dim_count = rect_dim_count + 1;
                    end
                end

                if rect_dim_count == size(g, 2)
                    bl = true;
                end
            end
        end

        function set_int_dx_(obj)
            g = obj.ss_grid_;
            if obj.ss_grid_is_rect_
                dx = zeros(1, size(g, 2));
                for lv = 1:size(g, 2)
                    rg = reshape(g(:, lv), sqrt(size(g(:, lv), 1)) * ones(1, 2));
                    drg1 = diff(rg, 1, 1);
                    drg2 = diff(rg, 1, 2);

                    dx(lv) = max(drg1(1), drg2(1));
                end

                %% For trapezoid rule, we penalize 1/2 per dimension at the endpoints
                % Endpoint iff one of the dimensions is xmin or xmax
                % Use max to implement OR and sum to count how many active dimensions
                trap_scale_factor = 2.^(-1 * ...
                    sum(obj.ss_grid_ == obj.xmin_ | ...
                        obj.ss_grid_ == obj.xmax_, 2));
                obj.int_dx_ = prod(dx) * trap_scale_factor;
            end
        end
    end

end