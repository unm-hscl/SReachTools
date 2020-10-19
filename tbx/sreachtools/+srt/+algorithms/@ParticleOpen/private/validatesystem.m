function validatesystem(obj, system)
% VALIDATESYSTEM Validate system.

arguments
    obj
    system (1, 1) srt.systems.StochasticSystem {mustBeNonempty}
end

supportedSystems = {'srt.systems.LtiSystem', ...
                    'srt.systems.LtvSystem'};

if ~ismember(class(system), supportedSystems)
  error('System is not supported.');
end

end
