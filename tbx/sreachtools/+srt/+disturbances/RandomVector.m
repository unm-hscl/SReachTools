classdef RandomVector
    properties (Access = protected)
        sample_fun_
    end

    methods
        function obj = RandomVector(varargin)
            p = inputParser();
            addOptional(p, 'sample_fun', @(n) [], ...
                @(x) isa(x, {'function_handle'}));
            parse(p, varargin{:});

            obj.sample_fun_ = p.Results.sample_fun;
        end

        function s = sample(obj, varargin)
            p = inputParser();
            addOptional(p, 'n', 1, @(x) validateattributes(x, {'numeric'}, ...
                {'scalar'}));
            parse(p, varargin{:});

            s = obj.sample_fun_(p.Results.n);
        end

        function rv = concat(obj)
            rv = srt.disturbances.RandomVector();
        end

        function rv = plus(A, B)
            if isa(A, 'srt.disturbances.RandomVector')
                rv = srt.disturbances.RandomVector(@(n) A.sample_fun_(n) + B);
            else
                rv = srt.disturbances.RandomVector(@(n) B.sample_fun_(n) + A);
            end
        end

        function rv = minus(A, B)
            if isa(A, 'srt.disturbances.RandomVector')
                rv = srt.disturbances.RandomVector(@(n) A.sample_fun_(n) - B);
            else
                rv = srt.disturbances.RandomVector(@(n) B.sample_fun_(n) - A);
            end
        end

        function rv = times(A, B)
            if isa(A, 'srt.disturbances.RandomVector')
                rv = srt.disturbances.RandomVector(@(n) B * A.sample_fun_(n));
            else
                rv = srt.disturbances.RandomVector(@(n) A * B.sample_fun_(n));
            end
        end

        function rv = mtimes(A, B)
            s = obj.sample();
            if isa(A, 'srt.disturbances.RandomVector')
                error(['Right matrix multiplication not supported for ', ...
                    'srt.disturbances.RandomVector']);
            else
                if size(A, 2) ~= length(s);
                    error(['Incorrect dimensions for matrix multiplication. ', ...
                        'Check that the number of columns in the matrix ', ...
                        'matches the lenght of a sample. To perform ', ...
                        'elementwise multiplication, use ''.*''.']);
                end

                rv = srt.disturbances.RandomVector(@(n) A * B.sample_fun_(n));
            end
        end

        function prob = getProbPolyhedron(obj, test_polyhedron, varargin)
        % Compute the probability of a random vector lying in a polyhedron
        % ====================================================================
        % The probability computation is done via Monte-Carlo (quasi Monte-Carlo
        % in Gaussian case), which is a guaranteed underapproximation (with a
        % probabilistic guarantee). The guarantee in the general case comes
        % via Hoeffding's inequality, and, in the Gaussian case, via the
        % confidence interval.
        %
        % General case: A distribution-independent lower bound on the number
        % of particles needed is obtained via Hoeffding's inequality.
        %
        % For \beta = exp(-2 * n_particles * \delta^2), Hoeffding's inequality
        % states that
        %
        %               Prob{X - E[X] \geq \delta} \leq \beta,
        %
        % where, X is the empirical average of a collection of n_particles
        % random variables bounded between [0,1], E[X] is the true mean, \delta
        % is the desired_accuracy, \beta is the failure risk (the probability of
        % the statement "X - E[X] \geq \delta" fails.
        %
        % In this function, we set the bounded random variables to be Bernoulli
        % random variables that is 1 when a realization of the random vector
        % is within the given test_polyhedron, and 0 otherwise.
        % Consequently, E[X] simplifies to the probability of the
        % realization to lie in the polyhedron.
        %
        % Given a failure risk \beta and desired accuracy \delta, we backcompute
        % the number of particles as the following
        %
        %               n_particles = -ln(\beta) / (2 * delta^2)
        %
        % This enforces the condition that empirical probability estimate X does
        % not overapproximate the true probability E[X] by more than delta.
        %
        % Gaussian case: We use Genz's algorithm to estimate the
        % probability. We increase the number of particles in the powers of
        % 10, till error estimate is within the desired tolerance.
        %
        % In both of the cases, to remain conservative in our estimate, we
        % subtract the desired_accuracy from the emperical probability estimate.
        % In other words, we use the lower bound of the confidence interval.
        % Thus, the probability estimate is guaranteed to be an
        % underapproximation.
        %
        % Inputs:
        % -------
        %   obj             - RandomVector object
        %   test_polyhedron - Polyhedron object (polytope whose probability of
        %                     occurrence is of interest)
        %   desired_accuracy- [Optional] Maximum absolute deviation from the
        %                     true probability estimate [Default: 1e-2]
        %
        % Outputs:
        % --------
        %   covar           - Probability of the random vector lying in the
        %                     given polytope
        % Notes:
        % ------
        % * Due to the inverse-square dependence on the desired_accuracy, we
        %   impose a hard lower bound of 1e-2 on the desired_accuracy. This
        %   leads to the requirement of 2e5 particles.
        % * We set the default failure risk as 2e-15.
        % * Ill-formed (mean/cov has Inf/NaN) Gaussian random vectors return
        %   zero probability.
        % * We compute the floor of the probability value based on the
        %   desired_accuracy to obtain a preferable lower bound on the
        %   probability.
        %
        % ====================================================================
        %
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        %
        %

            failure_risk = 2e-15;
            gauss_scale_factor = 10;
            if ~isempty(varargin)
                desired_accuracy = varargin{1};
                validateattributes(desired_accuracy, {'numeric'},...
                    {'scalar','>=', 1e-2, '<=', 1});
                if length(varargin) > 1
                    throwAsCaller(SrtInvalidArgsError('Too many inputs'));
                end
            else
                desired_accuracy = 1e-2;
            end

            % Check if we got a MPT's Polyhedron object of correct dimensions
            validateattributes(test_polyhedron, {'Polyhedron'}, {'nonempty'},...
                'RandomVector/getProbPolyhedron','test_polyhedron');
            if test_polyhedron.Dim ~= obj.dim
                throwAsCaller(SrtInvalidArgsError(['Mismatch in polytope ',...
                    'dimensions and random vector dimensions']));
            end

            switch obj.type
                case 'Gaussian'
                    if any(isinf(abs(obj.mean()))) || ...
                            any(isnan(abs(obj.mean()))) || ...
                            any(any(isinf(abs(obj.cov())))) || ...
                            any(any(isnan(abs(obj.cov()))))
                        % Return zero probability if ill-formed random vector
                        temp_probability = 0;
                    else
                        % Construct the half-space representation for qscmvnv
                        qscmvnv_lb = repmat(-Inf, ...
                            [size(test_polyhedron.A, 1),1]);
                        qscmvnv_coeff_matrix = test_polyhedron.A;
                        % We have the polytope given by Ax <= b, but qscmvnv
                        % expects a zero mean Gaussian. Hence, define x->x+mean
                        qscmvnv_ub = test_polyhedron.b - ...
                            test_polyhedron.A*obj.mean();
                        % Set up for an iterative call of Genz's particles;
                        % Increase particles till the error estimate is within
                        % threshold
                        n_particles = 10;
                        est_3sig_CI = 1;
                        while est_3sig_CI > desired_accuracy
                            try
                                % Call Genz's algorithm;
                                [temp_probability, est_3sig_CI] = qscmvnv( ...
                                    n_particles,obj.cov(), qscmvnv_lb, ...
                                    qscmvnv_coeff_matrix, qscmvnv_ub);
                            catch ME
                                disp(ME);
                                throw(SrtDevError(['Error in qscmvnv (', ...
                                    'Quadrature of multivariate Gaussian)']));
                            end
                            n_particles = n_particles * gauss_scale_factor;
                        end
                    end
                case 'UserDefined'
                    % By Hoeffding's inequality
                    n_particles = ceil(-log(failure_risk)/ ...
                        (2*desired_accuracy^2));

                    if n_particles > 2e5
                        warning('SReachTools:runtime', sprintf(['Number of ',...
                                'particles required: %1.2e. Consider ', ...
                                'increasing (relaxing) desired_accuracy!'], ...
                                n_particles));
                    end

                    mcarlo_sims = obj.getRealizations(n_particles);
                    count_contains = test_polyhedron.contains(mcarlo_sims);
                    temp_probability = sum(count_contains)/n_particles;
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['Probability computation not available for ',...
                         '%s-type random vectors'], obj.type)));
            end
            % Subtract the desired_accuracy to remain an underapproximation
            temp_probability = temp_probability - desired_accuracy;

            if temp_probability > 0
                % Rounding DOWN the integral to the desired accuracy
                prob = floor(temp_probability/desired_accuracy) * ...
                    desired_accuracy;
            else
                prob = 0;
            end
        end
    end

    % methods (Access = protected)
    %     function validatesample(obj)
    %         n = randi([1, 10]);
    %         s = obj.sample(n);
    %         if size(s, 2) ~= n
    %             error(['Sample function is not valid. Sampling functions ', ...
    %                 'must return vectors (nx1).']);
    %         end

    %         obj.n_ = size(s, 1);
    %     end
    % end
end
