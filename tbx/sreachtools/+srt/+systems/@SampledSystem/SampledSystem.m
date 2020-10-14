classdef SampledSystem < srt.systems.StochasticSystem
% SAMPLEDSYSTEM Samples from a stochastic system.

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

    end

    properties (Dependent)
        % X State samples.
        %
        %   Samples should be in column format, such that the matrix has
        %   dimensions [nxM], where n is the dimensionality of the state and M
        %   is the number of samples:
        %
        %   [ | | |     ]
        %   [ x x x ... ]
        %   [ | | |     ]
        X

        % U Input samples.
        %
        %   Samples should be in column format, such that the matrix has
        %   dimensions [mxM], where m is the dimensionality of the input and M
        %   is the number of samples:
        %
        %   [ | | |     ]
        %   [ u u u ... ]
        %   [ | | |     ]
        U

        % W Disturbance samples.
        %
        %   Samples should be in column format, such that the matrix has
        %   dimensions [pxM], where p is the dimensionality of the disturbance
        %   and M is the number of samples:
        %
        %   [ | | |     ]
        %   [ w w w ... ]
        %   [ | | |     ]
        W

        % Y Output samples.
        %
        %   Samples should be in column format, such that the matrix has
        %   dimensions [nxM], where n is the dimensionality of the state and M
        %   is the number of samples:
        %
        %   [ | | |     ]
        %   [ y y y ... ]
        %   [ | | |     ]
        Y
    end

    methods
        function obj = SampledSystem(varargin)
            % SAMPLEDSYSTEM Constructs an instance of the system.

            p = inputParser;
            addParameter(p, 'X', []);
            addParameter(p, 'U', []);
            addParameter(p, 'W', []);
            addParameter(p, 'Y', []);
            parse(p, varargin{:});

            M = size(p.Results.X, 2);

            % validateattributes(p.Results.U, {'numeric'}, {'ncols', M});
            % validateattributes(p.Results.W, {'numeric'}, {'ncols', M});
            % validateattributes(p.Results.Y, {'numeric'}, {'ncols', M});

            obj.num_samples_ = M;

            obj.x_samples_ = p.Results.X;
            obj.u_samples_ = p.Results.U;
            obj.w_samples_ = p.Results.W;
            obj.y_samples_ = p.Results.Y;

            obj.n_ = size(obj.x_samples_, 1);
            obj.m_ = size(obj.u_samples_, 1);
            obj.p_ = size(obj.w_samples_, 1);

            assert(obj.n_ == size(obj.y_samples_, 1));

        end
    end

    methods
        function len = length(obj)
            len = obj.num_samples_;
        end

        function X = get.X(obj)
            X = obj.x_samples_;
        end

        function U = get.U(obj)
            U = obj.u_samples_;
        end

        function W = get.W(obj)
            W = obj.w_samples_;
        end

        function Y = get.Y(obj)
            Y = obj.y_samples_;
        end
    end

end
