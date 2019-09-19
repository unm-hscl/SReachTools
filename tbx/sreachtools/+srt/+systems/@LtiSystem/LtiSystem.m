classdef LtiSystem < srt.systems.LtvSystem
    methods
        function obj = LtiSystem(varargin)
            p = inputParser();
            addOptional(p, 'A', [], @(x) validateattributes(x, ...
                {'numeric'}, {'square'}));
            addOptional(p, 'B', [], @(x) isa(x, 'numeric'));
            addOptional(p, 'F', [], ...
                @(x) srt.systems.LtiSystem.validateDisturbanceMatrix(x));
            addOptional(p, 'w', srt.disturbances.Empty(), ...
                @(x) isa(x, 'srt.disturbances.RandomVector'));
            parse(p, varargin{:});


            obj@srt.systems.LtvSystem(p.Results.A, p.Results.B, ...
                p.Results.F, p.Results.w)
        end

        function val = StateMatrix(obj, ~)
            val = obj.A;
        end

        function val = InputMatrix(obj, ~)
            val = obj.B;
        end

        function val = DisturbanceMatrix(obj, ~)
            val = obj.F;
        end

        function val = A(obj, ~)
            val = obj.A_;
        end

        function val = B(obj, ~)
            val = obj.B_;
        end

        function val = F(obj, ~)
            val = obj.F_;
        end
    end

    methods (Static)
        function valid = validateDisturbanceMatrix(F)
            if (isa(F, 'numeric') && length(size(F)) == 2) || ...
                (isa(F, 'char') && any(strcmp(F, {'InputMatrix', 'B'})))
                
                valid = true;
            else
                error(sprintf(['Expected input to be one of these types:\n\n', ...
                    'double, single, uint8, uint16, ', ...
                    'uint32, uint64, int8, int16, int32, int64\n\n', ...
                    'Or should a character array matching either:\n\n', ...
                    '''InputMatrix'', ''B''']));
            end
        end
    end
end