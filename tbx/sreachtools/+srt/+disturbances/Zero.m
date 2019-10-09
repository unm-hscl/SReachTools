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
            obj.sample_fun_ = @(n) zeros(obj.n_, n);
        end

        function rv = concat(obj)
            rv = srt.disturbances.Zero();
        end
    end
end