function validatesystem(sys, varargin)
% VALIDATESYSTEM Validate system object.

validateattributes(sys, {'srt.systems.StochasticSystem'}, {'nonempty'});
