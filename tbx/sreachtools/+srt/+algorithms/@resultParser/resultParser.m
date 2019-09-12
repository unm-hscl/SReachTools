classdef resultParser

    properties (Access = private)
        % Results struct.
        results_ struct
    end

    properties (Dependent)
        % RESULTS Result structure.
        Results
    end

    methods
        function results = get.Results(obj)
            % Return the results.
            results = obj.results_;
        end
    end

end
