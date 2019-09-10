classdef Empty < srt.spaces.Base
    methods
        function obj = Empty()
            obj@srt.spaces.Base()
        end
        
        function yn = contains(obj, x)
            yn = false;
        end

        function sp = concat(obj, time_horizon)
            sp = srt.spaces.Empty();
        end
    end
end