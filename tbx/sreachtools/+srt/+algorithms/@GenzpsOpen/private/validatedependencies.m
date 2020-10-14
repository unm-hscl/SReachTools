function validate_dependencies(obj)
% VALIDATE_DEPENDENCIES Checks if algorithm dependencies are present/valid.


v = ver;
%% Check for Global Optimization Toolbox (patternsearch)
has_patternsearch = any(strcmp(cellstr(char(v.Name)), ...
    'Global Optimization Toolbox'));
if ~has_patternsearch
    warning('SReachTools:setup',['''genzps-open'' option of ', ...
        'SReachPoint() function in SReachTools requires MATLAB''s ', ...
        'Global Optimization Toolbox.']);
end

end
