classdef Empty < srt.spaces.Base
    methods
        function obj = Empty()
            obj@srt.spaces.Base()
        end
        
        function yn = contains(obj, x)
            yn = false;
        end
    end
end