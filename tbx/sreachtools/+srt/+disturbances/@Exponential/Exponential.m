classdef Exponential < srt.disturbances.RandomVector
    properties (Dependent)
        % RATE Expoenential rate vector
        Rate

        % MEAN Exponential mean vector
        Mean
    end

    properties (Access = private)
        % RATE_ Exponential rate vector
        rate_

        % N_ Dimension of the disturbance vector
        n_
    end

    methods
        function obj = Exponential(rate)
            p = inputParser();
            addRequired(p, 'rate', @(x) validateattributes(x, {'numeric'}, ...
                {'positive', 'vector'}));
            parse(p, rate);

            obj@srt.disturbances.RandomVector();

            obj.rate_ = reshape(rate, [], 1);
            obj.n_ = length(obj.rate);

            obj.sample_fun = @(n) exprnd(repmat(obj.rate_, ...
                length(obj.rate_), n));
        end

        function val = get.Rate(obj)
            val = obj.rate_;
        end

        function val = get.Mean(obj)
            val = 1 ./ obj.Rate;
        end

        function y = pdf(obj, x)
            y = exppdf(x, obj.Rate);
        end

        function rv = concat(obj, time_horizon)
            rv = srt.disturbances.Exponential(kron(ones(time_horizon, 1), ...
                obj.Rate));
        end
    end
end