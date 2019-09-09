classdef (Abstract) Problem < handle
% PROBLEM Defines the interface for a problem.

  properties (Access = protected)
    time_horizon_(1, 1) double {mustBeInteger, mustBePositive} = 1

    % SAFE_SET_ Safe set must be a scalar.
    safe_set_(1, :) {mustBeValidSet} = Polyhedron()
    % TARGET_SET_ Target set must be a scalar.
    target_set_(1, 1) {mustBeValidSet} = Polyhedron()
  end

  properties (Dependent)
    % TIMEHORIZON The time horizon of the problem.
    TimeHorizon

    % SAFESET The safe set.
    %
    %   The SAFESET is the safe set of the system.
    %   Can be defined as either a Polyhedron, Tube, or as a Function.
    %
    %   See also: Polyhedron, Function
    SafeSet

    % TARGETSET The target set.
    %
    %   The TARGETSET is the target set of the system.
    %   Can be defined as either a Polyhedron, Tube, or as a Function.
    %
    %   See also: Polyhedron, Function
    TargetSet
  end

  methods
    function N = get.TimeHorizon(obj)
      N = obj.time_horizon_;
    end
    function set.TimeHorizon(N)
      obj.time_horizon_ = N;
    end

    function K = get.SafeSet(obj)
      K = obj.safe_set_;
    end
    function set.SafeSet(obj, K)
      if isa(K, function_handle)
        K = Function(K);
      end
      obj.safe_set_ = K;
    end

    function T = get.TargetSet(obj)
      T = obj.target_set_;
    end
    function set.TargetSet(obj, T)
      if isa(K, function_handle)
        K = Function(K);
      end
      obj.target_set_ = T;
    end
  end

  methods
    function tf = in_safe_set(obj, varargin)
      if isa(obj.safe_set_, 'Function')
        tf = obj.safe_set_.feval(varargin{:});
      else
        tf = obj.safe_set_.contains(varargin{:});
      end
    end
    function tf = in_target_set(obj, varargin)
      if isa(obj.target_set_, 'Function')
        tf = obj.target_set_.feval(varargin{:});
      else
        tf = obj.target_set_.contains(varargin{:});
      end
    end
  end

  methods (Abstract)
    tf = contains(obj, varargin)
  end

end

function mustBeValidSet(set)
  validateattributes(set, {'Polyhedron', 'Function'}, {'nonempty'});
end
