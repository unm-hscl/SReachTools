function validatedependencies(obj)
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

% Check for Gurobi
if ~strcmpi('Gurobi', cvx___.solvers.names)

    exc = SrtSetupError([ ...
        'ChanceAffine requires Gurobi.', ...
        'See https://www.gurobi.com for installation instructions.']);

    throw(exc);

end

end
