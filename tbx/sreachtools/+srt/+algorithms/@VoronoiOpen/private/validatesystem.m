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

% Must not be random vector.
% Must have Gaussian disturbance.

% 1. Get the correct options
% 2. Check if the system is stochastic
% 3. Check if the random vector is Gaussian
% 4. Check if the disturbance matrix is identity (Theory requirement,
%    else this formulation is not practical as it is stronger than
%    state feedback control law)

% if ~isa(sys.dist,'RandomVector')
%     throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
% end
% if ~strcmpi(sys.dist.type,'Gaussian')
%     throw(SrtInvalidArgsError(['SReachPointCcA requires Gaussian-',...
%         'perturbed LTV/LTI system']));
% end
% % Ensure that the disturbance matrix is identity
% err_string = ['Disturbance matrix must be an identity matrix.\nIf you ',...
%     'must have the disturbance matrix, redefine the disturbance random ',...
%     'vector.'];
% if sys.islti() && ~isequal(sys.dist_mat, eye(sys.dist_dim)) % Is F_k=I_k?
%     throwAsCaller(SrtInvalidArgsError(err_string));
% elseif sys.isltv()
%     not_equal_to_identity_flag = 0;
%     for t=0:time_horizon
%         if ~isequal(sys.dist_mat(t), eye(sys.dist_dim))
%             not_equal_to_identity_flag = 1;
%             break;
%         end
%     end
%     if not_equal_to_identity_flag
%         throwAsCaller(SrtInvalidArgsError(err_string));
%     end
% end

end
