classdef (Abstract) Algorithm < handle
% ALGORITHM Defines the interface for an algorithm.
%
%   Copyright 2019 Adam Thorpe

    properties (Access = private)
        verbose_(1, 1) double {mustBeNonnegative} = 0;
    end

    methods
        function obj = Algorithm(varargin)
            % ALGORITHM Construct an instance of the algorithm.

            p = inputParser;
            addParameter(p, 'verbose', 0);
            parse(p, varargin{:});

            obj.verbose_ = p.Results.verbose;

            % Validate algorithm dependencies.
            obj.validatedependencies();
        end
    end

    methods (Hidden, Access = protected)
        function validatedependencies(obj)
            % VALIDATEDEPENDENCIES Validate dependencies.
        end

        function validateproblem(obj, prob)
            % VALIDATEPROBLEM Validate problem.
        end

        function validatesystem(obj, sys)
            % VALIDATESYSTEM Validate system.
        end
    end

    methods
        function results = compute_point(obj, varargin)
            % COMPUTE_POINT Template method for algorithm.
            warning('compute_point not implemented for %s.', class(obj));
        end

        function results = compute_set(obj, varargin)
            % COMPUTE_SET Template method for algorithm.
            warning('compute_set not implemented for %s.', class(obj));
        end

        function results = compute_fwd(obj, varargin)
            % COMPUTE_FWD Template method for algorithm.
            warning('compute_fwd not implemented for %s.', class(obj));
        end
    end

    methods (Access = protected)
        function print_verbose(obj, level, varargin)
            % print_verbose Prints verbose messages to the command window.
            if obj.verbose_ >= level
                fprintf(varargin{:});
            end
        end
    end

end
