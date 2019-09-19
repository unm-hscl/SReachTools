classdef NdIntegrator < srt.systems.LtiSystem
    properties (Dependent)
        % SAMPLINGTIME Sampling time of the discrete integrator system
        SamplingTime
    end

    properties (Access = private)
        % T_ Sampling time of the discrete integrator system
        T_
    end

    methods
        function obj = NdIntegrator(n, T, varargin)
            p = inputParser();
            addRequired(p, 'n', @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
            addRequired(p, 'T', @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'positive'}));
            addOptional(p, 'F', [], ...
                @(x) srt.systems.LtiSystem.validateDisturbanceMatrix(x));
            addOptional(p, 'w', srt.disturbances.Empty(), ...
                @(x) isa(x, 'srt.disturbances.RandomVector'));
        
            parse(p, n, T, varargin{:});
        
            F = p.Results.F;
            w = p.Results.w;
        
            % anonymous function for getting the necessary matrix internals
            facT = @(t, n) t^n / factorial(n);
        
            % initialization
            A = eye(n);
            B = zeros(n, 1);
            
            % Populate the upper triangle of A and the entries of B
            for i = 1:n
                B(i) = facT(T, n-i+1);
                for j = i+1:n
                    A(i, j) = facT(T, j-i);
                end
            end

            obj = obj@srt.systems.LtiSystem(A, B, F, w);
        
            obj.T_ = T;
        end

        function val = get.SamplingTime(obj)
            val = obj.T_;
        end
    end
end
    