classdef DynamicProgramming < Algorithm
% DYNAMICPROGRAMMING Dynamic programming implementation.

    properties
        % State space grid.
        grid_(1, :) double {mustBeNumeric, mustBePositive}
        % Probability values.
        prob_ double {mustBeNumeric}
        % Probability levels.
        levels_ double {mustBeNumeric}
    end

    properties (Dependent)
        % GRID State space grid.
        Grid

    end

    methods
        function obj = DynamicProgramming(varargin)
            % DYNAMICPROGRAMMING Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@Algorithm(varargin{:});

        end
    end

    methods
        function grid = get.Grid(obj)

        end
        function set.Grid(obj, grid)

        end
    end

end
