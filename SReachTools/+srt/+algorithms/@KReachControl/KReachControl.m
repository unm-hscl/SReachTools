classdef KReachControl < Algorithm
% KReachControl Kernel distribution embeddings.

  properties (Access = private)
    % N State space dimensionality.
    n_(1, 1) double {mustBeNumeric}
    % M Input space dimensionality.
    m_(1, 1) double {mustBeNumeric}
    % P Disturbance space dimensionality.
    p_(1, 1) double {mustBeNumeric}

    % NUM_SAMPLES_ Total number of samples.
    num_samples_(1, 1) double {mustBeNumeric, mustBePositive} = 1024

    % X_SAMPLES_ Vector of state samples.
    x_samples_ double {mustBeNumeric}
    % U_SAMPLES_ Vector of input samples.
    u_samples_ double {mustBeNumeric}
    % W_SAMPLES_ Vector of disturbance samples.
    w_samples_ double {mustBeNumeric}
    % Y_SAMPLES_ Vector of output samples.
    y_samples_ double {mustBeNumeric}

    % SIGMA_ Sigma parameter to Gaussian kernel.
    sigma_(1, 1) double {mustBeNumeric, mustBePositive} = 0.1
    % LAMBDA_ Regularization parameter.
    lambda_(1, 1) double {mustBeNumeric, mustBePositive} = 1

    % VALUE_FUNCTIONS_ Computed value functions.
    value_functions_ double {mustBeNumeric}
  end

  % Hidden, Dependent properties. These can be accessed, but are hidden.
  properties (Hidden, Dependent)
    % X State samples.
    %
    %   If you pass samples into KReachControl manually, make sure they are in column
    %   format and that the matrix has the dimensions [nxM], where n is the
    %   dimensionality of the state and M is the number of samples:
    %
    %   [ | | |     ]
    %   [ x x x ... ]
    %   [ | | |     ]
    x

    % U Input samples.
    %
    %   If you pass samples into KReachControl manually, make sure they are in column
    %   format and that the matrix has the dimensions [mxM], where m is the
    %   dimensionality of the input and M is the number of samples:
    %
    %   [ | | |     ]
    %   [ u u u ... ]
    %   [ | | |     ]
    u

    % W Disturbance samples.
    %
    %   If you pass samples into KReachControl manually, make sure they are in column
    %   format and that the matrix has the dimensions [pxM], where p is the
    %   dimensionality of the disturbance and M is the number of samples:
    %
    %   [ | | |     ]
    %   [ w w w ... ]
    %   [ | | |     ]
    w

    % Y Output samples.
    %
    %   If you pass samples into KReachControl manually, make sure they are in column
    %   format and that the matrix has the dimensions [nxM], where n is the
    %   dimensionality of the state and M is the number of samples:
    %
    %   [ | | |     ]
    %   [ y y y ... ]
    %   [ | | |     ]
    y
  end

  properties (Dependent)
    % SIGMA Gaussian kernel bandwidth parameter.
    Sigma
    % LAMBDA Regularization parameter.
    Lambda
  end

  methods
    function obj = KReachControl(varargin)

      p = inputParser;
      if isa(varargin{1}, 'SampleGenerator');
        addRequired(p, 'gen');
      else
        addParameter(p, 'X', []);
        addParameter(p, 'U', []);
        addParameter(p, 'W', []);
        addParameter(p, 'Y', []);
      end
      addParameter(p, 'sigma', 0.1);
      addParameter(p, 'lambda', 1);
      parse(p, varargin{:});

      if isa(varargin{1}, 'SampleGenerator')
        gen = p.Results.gen;
        obj.x_samples_ = gen.generate_x_samples();
        obj.u_samples_ = gen.generate_u_samples();
        obj.w_samples_ = gen.generate_w_samples();
        obj.y_samples_ = gen.generate_y_samples();
      else
        obj.x_samples_ = p.Results.X;
        obj.u_samples_ = p.Results.U;
        obj.w_samples_ = p.Results.W;
        obj.y_samples_ = p.Results.Y;
      end

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
      cxx = obj.compute_norm(x);
      cxx = exp(-cxx/(2*sigma^2));
    end
    function cxy = compute_cross_covariance(x, y, sigma)
      % COMPUTE_CROSS_COVARIANCE Compute cross-covariance matrix.
      cxy = obj.compute_norm_cross(x, y);
      cxy = exp(-cxy/(2*sigma^2));
    end
  end

  methods
    function set.Sigma(obj, sigma)
      validateattributes(sigma, {'double'}, {'positive', 'scalar'});
      obj.sigma_ = sigma;
      obj.flag_param_ = true;
    end
    function set.Lambda(obj, lambda)
      validateattributes(lambda, {'double'}, {'positive', 'scalar'});
      obj.lambda_ = lambda;
      obj.flag_param_ = true;
    end
  end

end
