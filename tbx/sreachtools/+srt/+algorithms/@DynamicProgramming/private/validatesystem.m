function validatesystem(obj, sys)
% VALIDATESYSTEM Validate system.

% Cannot have large systems by default.

validateattributes(sys, {'srt.systems.StochasticSystem'}, {'nonempty'});

end
