classdef Exponential < srt.disturbances.RandomVector
    properties (SetAccess = private)
        rate
    end

    properties (Access = private)
        n_
    end

    methods
        function obj = Exponential(rate)
            obj.rate = rate;
            obj.n_ = length(obj.rate);
        end

        function s = sample(obj)
            s = exprnd(obj.rate);
        end

        function y = pdf(obj, x)
             y = mvnpdf(x, obj.rate);
        end

        function rv = concat(obj, time_horizon)
            rv = srt.disturbances.Exponential(kron(ones(time_horizon, 1), ...
                obj.rate));
        end
    end
end