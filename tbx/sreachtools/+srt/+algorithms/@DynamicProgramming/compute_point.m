function result = compute_point(obj, problem, sys, x0, varargin)

obj.validate_dependencies();
p = inputParser;
addRequired(p, 'problem', @obj.validate_problem);
addRequired(p, 'sys', @obj.validate_system);
addRequired(p, 'x0');
parse(p, problem, sys, x0, varargin{:});

end
