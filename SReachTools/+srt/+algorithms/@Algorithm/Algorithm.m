classdef (Abstract) Algorithm < handle
% ALGORITHM Defines the interface for an algorithm.

  methods (Abstract)
    validate_dependencies(obj)
    validate_problem(obj, problem)
    validate_system(obj, sys)
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

  methods (Sealed)
    function validate_arguments(varargin)
      % assert(problem.SafeSet.Dim == system.Dim);
      % assert(problem.TargetSet.Dim == system.Dim);
    end
  end

end
