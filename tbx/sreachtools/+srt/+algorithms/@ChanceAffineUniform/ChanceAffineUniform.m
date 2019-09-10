classdef ChanceAffineUniform < Algorithm
% CHANCEAFFINEUNIFORM Chance-constrained affine uniform.

    methods
        function obj = ChanceAffineUniform(varargin)
            % CHANCEAFFINEUNIFORM Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@Algorithm(varargin{:});

        end
    end
end
