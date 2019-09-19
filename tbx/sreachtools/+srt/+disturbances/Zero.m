classdef Zero < srt.disturbances.RandomVector
    properties (Access = private)
        n_
    end

    methods
        function obj = Zero(varargin)
            p = inputParser();
            addOptional(p, 'n', 1, @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer',  'positive'}));
            parse(p, varargin{:});

            obj@srt.disturbances.RandomVector()
            obj.n_ = p.Results.n;
        end

        function s = sample(obj)
            s = zeros(obj.n_, 1);
        end

        function rv = concat(obj)
            rv = srt.disturbances.Zero();
        end
    end
end