function validatesystem(obj, system)
% VALIDATESYSTEM Validate system.

supportedSystems = {'srt.systems.LtiSystem', ...
                    'srt.systems.LtvSystem'};

validateattributes(system, supportedSystems, {'scalar', 'nonempty'});

end
