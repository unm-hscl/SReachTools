classdef (Sealed) KernelEmbeddings < srt.algorithms.Algorithm
% KERNELEMBEDDINGS Kernel distribution embeddings.
%
%   Copyright 2019 Adam Thorpe

    properties (Access = private)

        % SIGMA_ Sigma parameter to Gaussian kernel.
        sigma_(1, 1) double {mustBeNumeric, mustBePositive} = 0.1
        % LAMBDA_ Regularization parameter.
        lambda_(1, 1) double {mustBeNumeric, mustBePositive} = 1

        % VALUE_FUNCTIONS_ Computed value functions.
        value_functions_ double {mustBeNumeric}
    end

    properties (Dependent)
        % SIGMA Gaussian kernel bandwidth parameter.
        Sigma
        % LAMBDA Regularization parameter.
        Lambda
    end

    methods
        function obj = KernelEmbeddings(varargin)
            % KERNELEMBEDDINGS Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            p = inputParser;
            addParameter(p, 'sigma', 0.1);
            addParameter(p, 'lambda', 1);
            parse(p, varargin{:});

            obj.sigma_ = p.Results.sigma;
            obj.lambda_ = p.Results.lambda;

        end
    end

    % Static methods.
    methods (Static, Hidden)
        function n = compute_norm(x)
            % COMPUTE_NORM Compute the norm.
            M = size(x, 2);
            n = zeros(M);

            for k = 1:size(x, 1)
                n = n + (repmat(x(k, :), [M, 1]) - repmat(x(k, :)', [1, M])).^2;
            end
        end
        function n = compute_norm_cross(x, y)
            % COMPUTE_CROSS_NORM Compute the cross norm.
            M = size(x, 2);
            T = size(y, 2);

            n = zeros(M, T);

            for k = 1:size(x, 1)
                n = n + (repmat(y(k, :), [M, 1]) - repmat(x(k, :)', [1, T])).^2;
            end
        end

        function n = normalize_beta(b)
            % NORMALIZE_BETA Normalize beta to ensure values are in [0, 1]
            n = b./sum(abs(b), 1);
        end

        function cxx = compute_autocovariance(x, sigma)
            % COMPUTE_AUTOCOVARIANCECOMPUTE Compute autocovariance matrix.
            cxx = srt.algorithms.KernelEmbeddings.compute_norm(x);
            cxx = exp(-cxx/(2*sigma^2));
        end
        function cxy = compute_cross_covariance(x, y, sigma)
            % COMPUTE_CROSS_COVARIANCE Compute cross-covariance matrix.
            cxy = srt.algorithms.KernelEmbeddings.compute_norm_cross(x, y);
            cxy = exp(-cxy/(2*sigma^2));
        end
    end

    methods
        function set.Sigma(obj, sigma)
            validateattributes(sigma, {'double'}, {'positive', 'scalar'});
            obj.sigma_ = sigma;
            obj.flag_param_ = true;
        end
        function sigma = get.Sigma(obj)
            sigma = obj.sigma_;
        end
        function set.Lambda(obj, lambda)
            validateattributes(lambda, {'double'}, {'positive', 'scalar'});
            obj.lambda_ = lambda;
            obj.flag_param_ = true;
        end
        function lambda = get.Lambda(obj)
            lambda = obj.lambda_;
        end
    end

end
