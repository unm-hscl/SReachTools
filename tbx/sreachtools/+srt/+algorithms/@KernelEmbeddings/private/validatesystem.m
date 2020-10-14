function validatesystem(obj, system)
% VALIDATESYSTEM Checks if the system is valid for the algorithm.

arguments
    obj
    system (1, 1) srt.systems.StochasticSystem {mustBeNonempty}
end

supportedSystems = {'srt.systems.SampledSystem'};

if ~ismember(class(system), supportedSystems)
  error('System is not supported.');
end

end
