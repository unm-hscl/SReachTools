classdef LtiSystem < srt.systems.LtvSystem
    methods
        function obj = LtiSystem(A, B, F, w)
            obj@srt.systems.LtvSystem(A, B, F, w)
        end

        function val = StateMatrix(obj, ~)
            val = obj.A;
        end

        function val = InputMatrix(obj, ~)
            val = obj.B;
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
end