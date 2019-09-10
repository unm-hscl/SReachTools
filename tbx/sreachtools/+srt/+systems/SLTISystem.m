classdef SLTISystem < srt.systems.SLTVSystem
    methods
        function obj = SLTISystem(A, B, F, w, U)
            obj@srt.systems.SLTVSystem(A, B, F, w, U)
        end

        function val = A(obj, k)
            val = obj.A_;
        end

        function val = B(obj, k)
            val = obj.B_;
        end

        function val = F(obj, k)
            val = obj.F_;
        end

        function xp = onestep(obj, x, varargin)
            p = inputParser();
            addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
            addOptional(p, 'u', [], ...
                @(x) validateattributes(x, {'numeric'}, {'vector'}));
            addParameter(p, 'UseMeanValue', false, @(x) islogical(x));
            parse(p, x, varargin{:});

            u = p.Results.u;

            Ax = obj.A * x;
            if p.Results.UseMeanValue
                Fw = obj.F * obj.w.mean();
            else
                Fw = obj.F * obj.w.sample();
            end
            Bu = obj.B * u;

            if isempty(Ax) && isempty(Bu) && isempty(Fw)
                xp = [];
            else
                if isempty(Ax) Ax = 0; end
                if isempty(Bu) Bu = 0; end
                if isempty(Fw) Fw = 0; end
                xp = Ax + Bu + Fw;
            end
        end

        function xm = meanstep(obj, x, u)
            p = inputParser()
            addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
            addOptional(p, 'u', [], ...
                @(x) validateattributes(x, {'numeric'}, {'vector'}));
            parse(p, x, varargin{:});

            xm = obj.onestep(k, x, p.Results.u, 'UseMeanValue', true)
        end
    end
end