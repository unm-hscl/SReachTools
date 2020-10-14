classdef (Abstract) Algorithm < handle
% ALGORITHM Defines the interface for an algorithm.
%
%   Copyright 2019 Adam Thorpe

    properties (Access = private)
        verbose_ (1, 1) double {mustBeNonnegative} = 0;
    end

    methods
        function obj = Algorithm(varargin)
            % ALGORITHM Construct an instance of the algorithm.

            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Verbose', 0);
            parse(p, varargin{:});

            obj.verbose_ = p.Results.Verbose;

            % Validate algorithm dependencies.
            obj.validatedependencies();
        end
    end

    methods (Hidden, Access = private)
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

    methods (Sealed, Access = protected)
        function print_verbose(obj, level, varargin)
            % print_verbose Prints verbose messages to the command window.
            if obj.verbose_ >= level
                fprintf(varargin{:});
            end
        end

        function plot_verbose(obj, level, cb)
            % print_verbose Prints verbose messages to the command window.
            if ~isa(cb, 'function_handle')
                error('Must provide a function handle.');
            end

            if obj.verbose_ >= level
                f = figure;
                ax = axes(f);
                cb(ax);
            end
        end
    end

end
