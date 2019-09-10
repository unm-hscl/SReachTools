classdef NdClosedBox < srt.spaces.Base
    properties (SetAccess = immutable)
        n

        lb

        ub
    end
    methods
        function obj = NdClosedBox(lb, ub)
            obj@srt.spaces.Base();

            validateattributes(lb, {'numeric'}, {'vector', 'real'})
            validateattributes(ub, {'numeric'}, {'vector', 'real'})

            obj.n  = length(lb);

            if size(lb, 2) > 1
                obj.lb = lb';
            else
                obj.lb = lb;
            end

            if size(ub, 2) > 1
                obj.ub = ub';
            else
                obj.ub = ub;
            end
        end

        function bl = contains(obj, v)
            if prod(size(v)) ~= obj.n
                error('Vector not of appropriate size')
            end

            bl = all(v >= obj.lb) && all(v <= obj.ub);
        end

        function plot(obj, varargin)
            p = inputParser();

            addParameter(p, 'Color', 'r', @(x) true);

            parse(p, varargin{:});

            lb = obj.lb;
            ub = obj.ub;

            if obj.n > 3
                error('Cannot plot NdBoxed of n > 3');
            elseif obj.n < 3
                patch([lb(1) ub(1) ub(1) lb(1)], [lb(2) lb(2) ub(2) ub(2)], ...
                    p.Results.Color);
            else
                patch([lb(1) ub(1) ub(1) lb(1)], ...
                    [lb(2) lb(2) ub(2) ub(2)], ...
                    [lb(2) lb(2) lb(2) lb(2)], p.Results.Color);

                patch([lb(1) ub(1) ub(1) lb(1)], ...
                    [lb(2) lb(2) ub(2) ub(2)], ...
                    [ub(2) ub(2) ub(2) ub(2)], p.Results.Color);

                patch([lb(1) ub(1) ub(1) lb(1)], ...
                    [lb(2) lb(2) lb(2) lb(2)], ...
                    [lb(2) lb(2) ub(2) ub(2)], p.Results.Color);

                patch([lb(1) ub(1) ub(1) lb(1)], ...
                    [ub(2) ub(2) ub(2) ub(2)], ...
                    [lb(2) lb(2) ub(2) ub(2)], p.Results.Color);

                patch([lb(1) lb(1) lb(1) lb(1)], ...
                    [lb(2) ub(2) ub(2) lb(2)], ...
                    [lb(2) lb(2) ub(2) ub(2)], p.Results.Color);

                patch([ub(1) ub(1) ub(1) ub(1)], ...
                    [lb(2) ub(2) ub(2) lb(2)], ...
                    [lb(2) lb(2) ub(2) ub(2)], p.Results.Color);
            end
        end

        function yn = isclosed(obj)
            yn = true;
        end

        function sp = concat(obj, time_horizon)
            sp = srt.spaces.NdClosedBox(kron(ones(time_horizon, 1), lb), ...
                kron(ones(time_horizon, 1), ub));
        end
    end
end