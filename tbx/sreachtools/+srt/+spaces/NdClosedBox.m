classdef NdClosedBox < srt.spaces.Base
    properties (SetAccess = immutable)
        n

        lb

        ub
    end
    methods
        function obj = NdClosedBox(lb, ub)
            obj@srt.spaces.Base();
            obj.n  = length(lb);
            obj.lb = lb;
            obj.ub = ub;
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
    end
end