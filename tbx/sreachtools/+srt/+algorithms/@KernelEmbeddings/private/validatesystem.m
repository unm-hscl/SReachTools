function validatesystem(obj, system)
% VALIDATESYSTEM Checks if the system is valid for the algorithm.

supportedSystems = {'srt.systems.SampledSystem'};

validateattributes(system, supportedSystems, {'scalar', 'nonempty'});

end
