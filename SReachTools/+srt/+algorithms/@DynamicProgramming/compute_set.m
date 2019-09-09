function result = compute_set(obj, problem, sys, varargin)

obj.validate_dependencies();
p = inputParser;
addRequired(p, 'problem', @obj.validate_problem);
addRequired(p, 'sys', @obj.validate_system);
parse(p, problem, sys, varargin{:});

end
