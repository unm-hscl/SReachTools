classdef (Abstract) Algorithm < handle
% ALGORITHM Defines the interface for an algorithm.
%
%   Copyright 2019 Adam Thorpe

    methods
        function obj = Algorithm(varargin)
            % ALGORITHM Construct an instance of the algorithm.

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
        function result = compute_point(obj, varargin)
            % COMPUTE_POINT Template method for algorithm.
            warning('compute_point not implemented for %s.', class(obj));
        end

        function result = compute_set(obj, varargin)
            % COMPUTE_SET Template method for algorithm.
            warning('compute_set not implemented for %s.', class(obj));
        end

        function result = compute_fwd(obj, varargin)
            % COMPUTE_FWD Template method for algorithm.
            warning('compute_fwd not implemented for %s.', class(obj));
        end
    end

end
