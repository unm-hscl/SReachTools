classdef FixedTurnRateDubinsCar < srt.systems.LtvSystem
    properties (Dependent)
        SamplingTime

        TurningRateSequence

        HeadingVector

        DisturbanceType
    end

    properties (Access = private)
        T_

        turn_seq_

        dist_type_
    end

    methods
        function obj = FixedTurnRateDubinsCar(varargin)

            p = inputParser();
            addParameter(p, 'TurningRateSequence', [], ...
                @(x) validateattributes(x, {'numeric'}, {'vector'}));
            addParameter(p, 'SamplingTime', [], @(x) validateattributes(x, ...
                {'numeric'}, {'scalar', 'positive'}));
            addParameter(p, 'DisturbanceMatrix', [], ...
                @(x) validateattributes(x, {'numeric'}, {'nonnan'}));
            addParameter(p, 'Disturbance', srt.disturbances.Empty(), ...
                @(x) isa(x, 'srt.disturbances.RandomVector'));
            addParameter(p, 'DisturbanceType', 'Additive', ...
                @(x) validatestring(x, {'Additive', 'Velocity'}));

            if isempty(p.Results.TurningRateSequence) || ...
                isempty(p.Results.SamplingTime) || ...
                isempty(p.Results.DisturbanceType)
                error(['Parameters ''TurningRateSequence'',', ...
                    '''SamplingTime'', and ''DisturbanceType'' must be ', ...
                    'specified.'])
            end

            T = p.Results.SamplingTime;

            turn_seq = reshape(p.Results.TurningRateSequence, 1, []);
            heading_vector = turn_seq(1) + T * cumsum([0, turn_seq(2:end)]);

            B = @(k) T * [cos(heading_vector(k+1)), sin(heading_vector(k+1))]';

            if strcmp(p.Results.DisturbanceType, 'Vector')
                dist_matrix = B;
            elseif strcmp(p.Results.DisturbanceType, 'Additive')
                dist_matrix = p.Results.DisturbanceMatrix;
            end

            % Superclass call
            obj = obj@srt.systems.LtvSystem(eye(2), B, dist_matrix, ...
                p.Results.Disturbance);

            obj.T_ = T;
            obj.turn_seq_ = turn_seq;
            obj.dist_type_ = p.Results.DisturbanceType;
        end

        function val = get.SamplingTime(obj)
            val = obj.T_;
        end

        function val = get.TurningRateSequence(obj)
            val = obj.turn_seq_;
        end

        function val = get.DisturbanceType(obj)
            val = obj.dist_type_;
        end

        function val = get.HeadingVector(obj)
            val = obj.turn_seq_(1) + obj.T_ * cumsum([0, obj.turn_seq_(2:end)]);
        end
    end
end