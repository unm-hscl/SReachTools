classdef Gaussian < srt.disturbances.RandomVector
    properties (SetAccess = private)
        mu
        sigma
    end
    properties (Access = private)
        n_
    end

    methods
        function obj = Gaussian(mu, sigma)
            obj@srt.disturbances.RandomVector()
            obj.mu    = mu;
            obj.sigma = sigma;

            obj.n_ = length(obj.mu);
        end

        function s = sample(obj, varargin)
            p = inputParser();

            addOptional(p, 'n', 1, @(x) validateattributes(x, {'numeric'}, ...
                {'scalar'}));

            parse(p, varargin{:});

            s = mvnrnd(obj.mu, obj.sigma, p.Results.n)';
        end

        function v = pdf(obj, x)
            v = mvnpdf(x, obj.mu, obj.sigma);
        end

        function v = cov(obj)
            v = obj.sigma;
        end

        function ss = sampleSpace(obj)
            ss = srt.spaces.Rn(obj.n_);
        end

        function rv = concat(obj, time_horizon)
            rv = srt.disturbances.Gaussian( ...
                kron(ones(time_horizon, 1), obj.mu), ...
                kron(eye(time_horizon), obj.sigma));
        end
    end
end