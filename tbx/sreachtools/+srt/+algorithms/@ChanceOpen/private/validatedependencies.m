function validate_dependencies(obj)
% VALIDATE_DEPENDENCIES Checks if algorithm dependencies are present/valid.

% Requires CVX.

% Check for CVX.
global cvx___

if ~isfield(cvx___, 'loaded') || ~cvx___.loaded

    exc = SrtSetupError([ ...
        'ChanceAffine requires CVX.', ...
        'See http://cvxr.com/cvx/download/ for installation instructions.']);

    throw(exc);

end

end
