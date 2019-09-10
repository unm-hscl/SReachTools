classdef ChanceAffine < Algorithm
% CHANCEAFFINE Chance-constrained affine.

    methods
        function obj = ChanceAffine(varargin)
            % CHANCEAFFINE Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@Algorithm(varargin{:});
            
        end
    end

end
