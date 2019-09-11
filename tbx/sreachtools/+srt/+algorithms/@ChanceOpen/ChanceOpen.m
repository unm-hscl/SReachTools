classdef ChanceOpen < srt.algorithms.Algorithm
% CHANCEOPEN Chance-constrained open-loop.

    methods
        function obj = ChanceOpen(varargin)
            % CHANCEOPEN Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

        end
    end

end
