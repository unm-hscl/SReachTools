function validatesystem(obj, system)
% VALIDATESYSTEM Checks if the system is valid for the algorithm.

validateattributes(system, {'srt.systems.StochasticSystem'}, {'nonempty'});

supportedSystems = {'srt.systems.SampledSystem'};

if ~ismember(class(system), supportedSystems)
  error('System is not supported.');
end

end
