classdef Gaussian < srt.disturbances.RandomVector
    properties (Dependent)
        % Mu Gaussian mean vector
        Mu

        % SIGMA Gaussian covariance matrix
        Sigma

        % MEAN Gaussian mean vector
        Mean
    end

    properties (Access = private)
        % MU_ Gaussian mean vector
        mu_

        % SIGMA_ Gaussian covariance matrix
        sigma_

        % N_ Dimension of the random vector
        n_
    end

    methods
        function obj = Gaussian(mu, sigma)
            obj@srt.disturbances.RandomVector();

            p = inputParser();
            addRequired(p, 'mu', @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
            addRequired(p, 'sigma', @(x) validateattributes(x, {'numeric'}, ...
                {'square'}));
            parse(p, mu, sigma);

            obj.mu_ = reshape(mu, [], 1);
            obj.sigma_ = sigma;

            obj.n_ = length(obj.mu_);
            obj.sample_fun_ = @(n) mvnrnd(obj.mu_, obj.sigma_, n);
        end
    
        function val = get.Mu(obj)
            val = obj.mu_
        end

        function val = get.Sigma(obj)
            val = obj.sigma_
        end

        function val = get.Mean(obj)
            val = obj.Mu;
        end

        function s = sample(obj, varargin)
            p = inputParser();
            addOptional(p, 'n', 1, @(x) validateattributes(x, {'numeric'}, ...
                {'scalar'}));
            parse(p, varargin{:});

            s = mvnrnd(obj.mu_, obj.sigma_, p.Results.n)';
        end

        function v = pdf(obj, x)
            v = mvnpdf(x, obj.mu_, obj.sigma_);
        end

        function ss = sampleSpace(obj)
            ss = srt.spaces.Rn(obj.n_);
        end

        function rv = concat(obj, time_horizon)
            rv = srt.disturbances.Gaussian( ...
                kron(ones(time_horizon, 1), obj.mu_), ...
                kron(eye(time_horizon), obj.sigma_));
        end
    end
end