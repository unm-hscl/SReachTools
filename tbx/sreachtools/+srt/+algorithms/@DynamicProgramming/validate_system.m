function validate_system(obj, system)
% VALIDATE_SYSTEM Checks if the system is valid for the algorithm.

validateattributes(system, {'StochasticSystem'}, {'nonempty'});

supportedSystems = {'StochasticLTISystem'};

if ~ismember(class(system), supportedSystems)
  error('System is not supported.');
end

end
