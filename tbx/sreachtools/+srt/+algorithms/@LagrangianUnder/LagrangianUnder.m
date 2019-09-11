classdef LagrangianUnder < srt.algorithms.Algorithm
% LAGRANGIANUNDER Lagrangian underapproximation.

    methods
        function obj = LagrangianUnder(varargin)
            % LAGRANGIANUNDER Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

        end
    end

end
