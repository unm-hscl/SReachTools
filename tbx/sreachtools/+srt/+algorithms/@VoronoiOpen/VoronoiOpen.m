classdef VoronoiOpen < srt.algorithms.Algorithm
% VORONOIOPEN Voronoi partition control.

    methods
        function obj = VoronoiOpen(varargin)
            % VORONOIOPEN Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

        end
    end

end
