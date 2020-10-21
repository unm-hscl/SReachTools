classdef Rn < srt.spaces.NdOpenBox
    methods
        function obj = Rn(n)
            lb = -Inf * ones(n, 1);
            ub = Inf * ones(n, 1);
            obj@srt.spaces.NdOpenBox(lb, ub);   
        end

        function yn = contains(obj, v)
            yn = prod(size(v)) == obj.n && any(size(v) == 1) && ~any(isinf(v));
        end

        function sp = concat(obj, time_horizon)
            sp = srt.spaces.Rn(time_horizon * obj.n);
        end
    end
end