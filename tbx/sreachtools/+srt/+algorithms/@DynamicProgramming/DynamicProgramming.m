classdef DynamicProgramming < srt.algorithms.Algorithm
% DYNAMICPROGRAMMING Dynamic programming implementation.

    properties (Dependent)


        Probability
    end

    properties
        % State space grid.
        ss_grid_(1, :) double {mustBeNumeric}

        % Input space grid.
        is_grid_(1, :) double {mustBeNumeric}

        % Probability values at state space grid points.
        prob_ double {mustBeNumeric}
    end

    properties (Dependent)
        % STATESPACEGRID State space grid.
        StateSpaceGrid

        % INPUTSPACEGRID Input space grid.
        InputSpaceGrid
    end

    methods
        function obj = DynamicProgramming(varargin)
            % DYNAMICPROGRAMMING Construct an instance of the algorithm.

            p = inputParser();
            addParameter(p, 'StateSpaceGrid', []);
            addParameter(p, 'InputSpaceGrid', []);
            addParameter(p, 'verbose', 0, @(x) validateattributes(x, ...
                {'numeric'}, {'scalar', 'integer', 'nonnegative'}));

            parse(p, varargin{:});

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            if isempty(p.Results.StateSpaceGrid) || ...
                isempty(p.Results.InputSpaceGrid)

                error(['Must specifty StateSpaceGrid and InputSpaceGrid ', ...
                    'parameters.']);
            end

            obj.ss_grid_ = obj.check_and_set_grid_(p.Results.StateSpaceGrid);
            obj.is_grid_ = obj.check_and_set_grid_(p.Results.InputSpaceGrid);

        end
    end

    methods
        function grid = get.Grid(obj)

        end
        function set.Grid(obj, grid)

        end
    end

    methods (Access = private)
        function g = check_and_set_grid_(obj, grid)
            if isnumeric(grid) && ismatrix(grid)
                g = grid;
                return;
            end

            if isnumeric(grid) && isvector(grid)
                g = reshape(grid, [], 1);
                return;
            end

            if iscell(grid)
                % Three ways of specifying the grid using a cell array
                %   1) {integer=n, vector}
                %       Builds with [X1, ..., Xn] = ndgrid(vector)
                %   2) {vector1, vector2, ...}
                %       Builds with [X1, X2, ...] = ndgrid(vector1, vector2, ...)
                %   3) {X1, X2, ...}
                %       Assumes X1, X2, ... are formed from meshgrid or ndgrid.
                %       Just reshapes
                % 

                valid_grid_input = true;

                if isnumeric(grid{1}) && isscalar(grid{1})
                    if length(grid) == 2 && isnumeric(grid{2}) ...
                        && isvector(grid{2})
                        scalar_vector = true;

                    else
                        scalar_vector = false;
                    end
                end

                if ~scalar_vector
                    all_vectors = true;
                    all_grids   = true;

                    for sub = grid
                        if isnumeric(sub) && isvector(sub)
                            all_vectors = all_vectors * true;
                            all_grids   = all_grids * false;
                        elseif isnumeric(sub) && length(size(sub)) == length(grid)
                            all_vectors = all_vectors * false;
                            all_grids   = all_grids * true;
                        else
                            all_vectors = all_vectors * false;
                            all_grids   = all_grids * false;
                            break;
                        end
                    end
                end

                if scalar_vector
                    M = cell(1, grid{1});
                    [M{:}] = ndgrid(grid{2});
                elseif all_vectors
                    M = cell(1, grid{1});
                    [M{:}] = ndgrid(grid{:});
                elseif all_grids
                    M = grid;
                else
                    error('Invalid grid input: Please see help.');
                end

                g = zeros(numel(M{1}), length(M));
                for i = 1:length(M)
                    g(i, :) = reshape(M{i}, [], 1);
                end
                    
            end
        end
    end

end