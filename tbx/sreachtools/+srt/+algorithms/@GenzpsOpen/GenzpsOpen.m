classdef GenzpsOpen < srt.algorithms.Algorithm
% GENZPSOPEN Genzps open.

    properties (Access = private)

        desired_accuracy (1, 1) double {mustBeLessThan(desired_accuracy, 1), ...
            mustBeGreaterThanOrEqual(desired_accuracy, 1E-2)} = 5E-2

        PSoptions (1, 1) = psoptimset('display','off');

        thresh (1, 1) double {mustBePositive, ...
            mustBeLessThanOrEqual(thresh, 1)} =  1

    end

    methods
        function obj = GenzpsOpen(varargin)
            % GENZPSOPEN Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            p = inputParser;
            p.KeepUnmatched = true;

            addParameter(p, 'desired_accuracy', 5E-2);
            addParameter(p, 'PSoptions', psoptimset('display','off'));
            addParameter(p, 'thresh', 1);

            parse(p, varargin{:});

            obj.desired_accuracy = p.Results.desired_accuracy;
            obj.PSoptions = p.Results.PSoptions;
            obj.thresh = p.Results.thresh;

        end
    end

end
