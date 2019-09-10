function validatesystem(obj, sys)
% VALIDATESYSTEM Validate system.

% Must not be random vector.
% Must have Gaussian disturbance.

% 1. Get the correct options
% 2. Check if the system is stochastic
% 3. Check if the random vector is Gaussian

% if ~isa(sys.dist,'RandomVector')
%     throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
% end
% if ~strcmpi(sys.dist.type,'Gaussian')
%     throw(SrtInvalidArgsError(['SReachPointCcO requires Gaussian-',...
%         'perturbed LTV/LTI system']));
% end

end
