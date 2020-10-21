classdef Empty < srt.disturbances.RandomVector
    methods
        function obj = Empty()
            obj@srt.disturbances.RandomVector()
        end

        function s = sample(obj)
            s = [];
        end

        function rv = concat(obj)
            rv = srt.disturbances.Empty();
        end
    end
end