classdef Zero < srt.disturbances.RandomVector
    methods
        function obj = Zero()
            obj@srt.disturbances.RandomVector()
        end

        function s = sample(obj)
            s = 0;
        end

        function rv = concat(obj)
            rv = srt.disturbances.Zero();
        end
    end
end