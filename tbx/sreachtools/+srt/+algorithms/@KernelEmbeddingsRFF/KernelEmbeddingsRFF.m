classdef (Sealed) KernelEmbeddingsRFF < srt.algorithms.Algorithm
% KERNELEMBEDDINGSRFF Kernel distribution embeddings (RFF).
%
%   Kernel distribution embeddings using random fourier features.
%
%   Copyright 2019 Adam Thorpe

    properties (Access = private)

        % Sigma parameter to Gaussian kernel.
        sigma_(1, 1) double {mustBeNumeric, mustBePositive} = 0.1
        % Regularization parameter.
        lambda_(1, 1) double {mustBeNumeric, mustBePositive} = 1

        % Computed value functions.
        value_functions_ double {mustBeNumeric}

        % Fourier dimension.
        D_
    end

    properties (Dependent)
        % SIGMA Gaussian kernel bandwidth parameter.
        Sigma
        % LAMBDA Regularization parameter.
        Lambda

        % D Fourier dimension.
        D
    end

    methods
        function obj = KernelEmbeddingsRFF(varargin)
            % KERNELEMBEDDINGSRFF Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'sigma', 0.1);
            addParameter(p, 'lambda', 1);
            addParameter(p, 'D', 1);
            parse(p, varargin{:});

            obj.sigma_ = p.Results.sigma;
            obj.lambda_ = p.Results.lambda;

            obj.D_ = p.Results.D;

        end
    end

    % Static methods.
    methods (Static, Hidden)
        function d = compute_diff(w, x)
            % COMPUTE_DIFF Compute the difference.
            D = size(w, 2);
            M = size(x, 2);

            d = zeros(D);

            wx = w.'*x;

            for k = 1:D
                d = d + repmat(wx(k, :), [D, 1]) - repmat(wx(k, :)', [1, D]);
            end
        end

        function n = compute_norm(x)
            % COMPUTE_NORM Compute the norm.
            M = size(x, 2);
            n = zeros(M);

            for k = 1:size(x, 1)
                n = n + (repmat(x(k, :), [M, 1]) - repmat(x(k, :)', [1, M]));
            end
        end
        function n = compute_norm_cross(x, y)
            % COMPUTE_CROSS_NORM Compute the cross norm.
            M = size(x, 2);
            T = size(y, 2);

            n = zeros(M, T);

            for k = 1:size(x, 1)
                n = n + (repmat(y(k, :), [M, 1]) - repmat(x(k, :)', [1, T]));
            end
        end

        function n = normalize_beta(b)
            % NORMALIZE_BETA Normalize beta to ensure values are in [0, 1]
            n = b./sum(abs(b), 1);
        end

        function cxx = compute_autocovariance(w, x, b)
            % COMPUTE_AUTOCOVARIANCECOMPUTE Compute autocovariance matrix.
            cxx = srt.algorithms.KernelEmbeddings.compute_norm(x);
            cxx = cos(w*cxx + b);
        end
        function cxy = compute_cross_covariance(w, x, y, b)
            % COMPUTE_CROSS_COVARIANCE Compute cross-covariance matrix.
            cxy = srt.algorithms.KernelEmbeddings.compute_norm_cross(x, y);
            cxy = cos(w*cxy + b);
        end
    end

    methods
        function set.Sigma(obj, sigma)
            validateattributes(sigma, {'double'}, {'positive', 'scalar'});
            obj.sigma_ = sigma;
        end
        function sigma = get.Sigma(obj)
            sigma = obj.sigma_;
        end
        function set.Lambda(obj, lambda)
            validateattributes(lambda, {'double'}, {'positive', 'scalar'});
            obj.lambda_ = lambda;
        end
        function lambda = get.Lambda(obj)
            lambda = obj.lambda_;
        end
        function set.D(obj, D)
            validateattributes(D, {'double'}, {'positive', 'scalar'});
            obj.D_ = D;
        end
        function D = get.D(obj)
            D = obj.D_;
        end
    end

end
