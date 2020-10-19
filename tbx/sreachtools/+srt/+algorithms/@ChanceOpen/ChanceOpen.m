classdef ChanceOpen < srt.algorithms.Algorithm
% CHANCEOPEN Chance-constrained open-loop.

    properties (Access = private)

        % Accuracy of piecewise-affine approximation of norminvcdf
        pwa_accuracy (1, 1) double {mustBePositive} = 1E-3

    end

    methods
        function obj = ChanceOpen(varargin)
            % CHANCEOPEN Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            p = inputParser;
            p.KeepUnmatched = true;

            addParameter(p, 'pwa_accuracy', 1E-3);

            parse(p, varargin{:});

            obj.pwa_accuracy = p.Results.pwa_accuracy;

        end
    end

end
