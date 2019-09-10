classdef Rn < srt.spaces.NdOpenBox
    methods
        function obj = Rn(n)
            lb = -Inf * ones(1, n);
            ub = Inf * ones(1, n);
            obj@srt.spaces.NdOpenBox(lb, ub);   
        end

        function yn = contains(obj, v)
            yn = prod(size(v)) == obj.n && any(size(v) == 1) && ~any(isinf(v));
        end
    end
end