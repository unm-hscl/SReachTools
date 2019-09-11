classdef (Abstract) Problem < handle
% PROBLEM Defines the interface for a problem.

    properties (Access = protected)
        time_horizon_(1, 1) double {mustBeInteger, mustBePositive} = 1

        % SAFE_TUBE_ Target set must be a scalar.
        safe_tube_(1, 1) {mustBeValidSet};
        % TARGET_TUBE_ Target set must be a scalar.
        target_tube_(1, 1) {mustBeValidSet};
    end

    properties (Dependent)
        % TIMEHORIZON The time horizon of the problem.
        TimeHorizon

        % SAFETUBE The target tube.
        %
        %   The TARGETTUBE is the target tube for the problem.
        %   Can be defined as either a Polyhedron, Function, or as a Tube.
        %
        %   See also: Polyhedron, Function, Tube
        SafeTube

        % TARGETTUBE The target tube.
        %
        %   The TARGETTUBE is the target tube for the problem.
        %   Can be defined as either a Polyhedron, Function, or as a Tube.
        %
        %   See also: Polyhedron, Function, Tube
        TargetTube
    end

    methods
        function N = get.TimeHorizon(obj)
            N = obj.time_horizon_;
        end
        function set.TimeHorizon(N)
            obj.time_horizon_ = N;
        end

        function T = get.SafeTube(obj)
            T = obj.safe_tube_;
        end
        function set.SafeTube(obj, T)
            if isa(K, function_handle)
                K = Function(K);
            end
            obj.safe_tube_ = T;
        end

        function T = get.TargetTube(obj)
            T = obj.target_tube_;
        end
        function set.TargetTube(obj, T)
            if isa(K, function_handle)
                K = Function(K);
            end
            obj.target_tube_ = T;
        end
    end

    methods
        function tf = in_target_tube(obj, varargin)
            if isa(obj.target_tube_, 'Function')
                tf = obj.target_tube_.feval(varargin{:});
            else
                tf = obj.target_tube_.contains(varargin{:});
            end
        end
    end

    methods (Abstract)
        tf = contains(obj, varargin)
    end

end

function mustBeValidSet(set)
    % isa(set, {'Polyhedron', 'Tube', 'Function'});
end
