function validate_dependencies(obj)
% VALIDATE_DEPENDENCIES Checks if algorithm dependencies are present/valid.

% Ensure that patternsearch is installed
v = ver;
has_patternsearch = any(strcmp(cellstr(char(v.Name)), ...
    'Global Optimization Toolbox'));
if ~has_patternsearch
    exc = SrtSetupError(['SReachPoint with ''genzps-open'' ', ...
        'option needs MATLAB''s Global Optimization Toolbox.']);
    throw(exc);
end

end
