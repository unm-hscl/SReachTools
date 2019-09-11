classdef LagrangianOver < srt.algorithms.Algorithm
% LAGRANGIANOVER Lagrangian overapproximation.

    methods
        function obj = LagrangianOver(varargin)
            % LAGRANGIANOVER Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

        end
    end

end
