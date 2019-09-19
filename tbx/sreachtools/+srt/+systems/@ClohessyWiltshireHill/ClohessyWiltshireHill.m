classdef ClohessyWiltshireHill < srt.systems.LtiSystem
    properties (Dependent)
        % SAMPLINGTIME Sampling time of the system
        SamplingTime

        % ORBITALRADIUS Deputy (satellite) orbital radius
        OrbitalRadius

        % CELESTIALBODYMASS Mass of the celestial body
        CelestialBodyMass

        % DEPUTYMASS Mass of the deputy (satellite) [kg]
        DeputyMass

        % MEANMOTION Mean motion of the deputy (satellite) [1/s]
        MeanMotion

        % MEANANOMALY Mean anomaly of the deputy (satellite)
        MeanAnomaly

        % GRAVITATIONALPARAMETER Gravitational parameter of the celesital body
        % [m^3 / s^2]
        GravitationalParameter
    end

    properties (Access = private)
        % T_ Sampline time of the system [s]
        T_
        
        % ORBITAL_RADIUS_ Deputy (satellite) orbital radius [km]
        orbital_radius_

        % CELESTIAL_MASS_ Mass of the celestial body [kg]
        celestial_mass_

        % DEPUTY_MASS_ Mass of the deputy (satellite) [kg]
        deputy_mass_
    end
    
    properties (Constant, Access = private)
        % G_ Universal gravitational constant [m^3 / (kg * s^2)]
        G_ = 6.673e-11;

        % EARTH_RADIUS_ Equitorial Radius of the earth [km]
        EARTH_RADIUS_ = 6378.1;

        % EARTH_MASS_ Mass of the Earth [kg]
        EARTH_MASS_ = 5.9472e24;
    end

    methods
        function obj = ClohessyWiltshireHill(varargin)
            import srt.systems.ClohessyWiltshireHill
            p = inputParser();
            addParameter(p, 'SamplingTime', [], @(x) validateattributes(x, ...
                {'numeric'}, {'scalar', 'positive', 'real'}));
            addParameter(p, 'OrbitalRadius', ...
                850 + ClohessyWiltshireHill.EARTH_RADIUS_, ...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'positive', 'real'}));
            addParameter(p, 'CelestialBodyMass', ...
                ClohessyWiltshireHill.EARTH_MASS_, ...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'positive', 'real'}));
            addParameter(p, 'DeputyMass', 300, @(x) validateattributes(x, ...
                {'numeric'}, {'scalar', 'positive', 'real'}));
            addParameter(p, 'Dimension', 6, @(x) x == 4 | x == 6);
            addParameter(p, 'DisturbanceMatrix', [], ...
                @(x) srt.systems.LtiSystem.validateDisturbanceMatrix(x));
            addParameter(p, 'Disturbance', srt.disturbances.Empty(), ...
                @(x) isa(x, 'srt.disturbances.RandomVector'));
            parse(p, varargin{:});

            if isempty(p.Results.SamplingTime)
                error('SamplingTime property must be specified');
            end

            grav_param = ClohessyWiltshireHill.G_ * ...
                p.Results.CelestialBodyMass / 1000^3;
            mean_motion = sqrt(grav_param / p.Results.OrbitalRadius^3);
            [A, B] = ClohessyWiltshireHill.getCwhStateAndInputMatrices( ...
                p.Results.SamplingTime, mean_motion, p.Results.DeputyMass);

            if p.Results.Dimension == 4
                % reduce matrices
                A = A([1:2, 4:5], [1:2, 4:5]);
                B = B([1:2, 4:5], 1:2);
            end
            
            obj = obj@srt.systems.LtiSystem(A, B, ...
                p.Results.DisturbanceMatrix, p.Results.Disturbance)

            obj.T_ = p.Results.SamplingTime;
            obj.orbital_radius_ = p.Results.OrbitalRadius;
            obj.celestial_mass_ = p.Results.CelestialBodyMass;
            obj.deputy_mass_ = p.Results.DeputyMass;
        end

        function val = get.SamplingTime(obj)
            val = obj.T_;
        end

        function val = get.OrbitalRadius(obj)
            val = obj.orbital_radius_;
        end

        function val = get.CelestialBodyMass(obj)
            val = obj.celestial_mass_;
        end

        function val = get.DeputyMass(obj)
            val = obj.deputy_mass_;
        end

        function val = get.GravitationalParameter(obj)
            val = obj.G_ * obj.celestial_mass_;
        end

        function val = get.MeanMotion(obj)
            val = sqrt(obj.GravitationalParameter / ...
                (obj.orbital_radius_ * 1000^3));
        end

        function val = get.MeanAnomaly(obj)
            val = obj.MeanMotion * obj.T_;
        end
    end

    methods (Static, Access = private)
        function [A, B] = getCwhStateAndInputMatrices(T, n, deputy_mass)
            
                % Continuous-time LTI CWH unforced dynamics e^{A_{cts}t}
            eAt = @(t) [ 
                4 - 3 * cos(n * t), ...
                    0, ...
                    0, ...
                    (1/n) * sin(n * t), ...
                    (2/n) * (1 - cos(n * t)), ...
                    0; ...
                6 * (sin(n * t) - n * t), ...
                    1, ...
                    0, ...
                    -(2/n) * (1 - cos(n * t)), ...
                    (1/n) * (4*sin(n * t) - 3*n * t), ...
                    0; ...
                0, ...
                    0, ...
                    cos(n * t), ...
                    0, ...
                    0, ...
                    (1/n) * sin(n * t); ...
                3 * n * sin(n * t), ...
                    0, ...
                    0, ...
                    cos(n * t), ...
                    2 * sin(n * t), ...
                    0; ...
                -6 * n * (1 - cos(n * t)), ...
                    0, ...
                    0, ...
                    -2 * sin(n * t), ...
                    4 * cos(n * t) - 3, ...
                    0; ...
                0, ...
                    0, ...
                    -n * sin(n * t), ...
                    0, ...
                    0, ...
                    cos(n * t);
            ];

            % Discrete-time state matrix is Phi(T_s) for sampling time T_s since the
            % system is time-invariant
            A = eAt(T);
        
            % Continuous-time input matrix B_{cts}
            B_cts = 1 / deputy_mass * [zeros(3); eye(3)];

            % Discrete-time input matrix is (\int_0^T e^{A_{cts}\tau} d\tau) B_cts
            B = integral(eAt, 0, T, 'ArrayValued', true) * B_cts;
        end
    end
end