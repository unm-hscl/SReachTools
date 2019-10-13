classdef (Abstract) Problem < handle
% PROBLEM Defines the interface for a problem.

    properties (Access = protected)
        constraint_tube_(1, 1) {mustBeValidTube};
        target_tube_(1, 1) {mustBeValidTube};
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
        ConstraintTube

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
            N = length(obj.constraint_tube_);
        end

        function T = get.ConstraintTube(obj)
            T = obj.constraint_tube_;
        end
        function set.ConstraintTube(obj, T)
            if isa(T, 'function_handle')
                T = Function(T);
            end
            obj.constraint_tube_ = T;
        end

        function T = get.TargetTube(obj)
            T = obj.target_tube_;
        end
        function set.TargetTube(obj, T)
            if isa(T, 'function_handle')
                T = Function(T);
            end
            obj.target_tube_ = T;
        end
    end

end

function mustBeValidTube(t)
% MUSTBEVALIDTUBE Verify tube is valid.

% Sets must be valid Tube.
% validateattributes(t, {'srt.Tube'}, {'scalar'});

end
