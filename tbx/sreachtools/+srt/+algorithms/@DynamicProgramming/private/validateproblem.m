function validateproblem(obj, problem)
% VALIDATEPROBLEM Checks if the problem is valid for the algorithm.

validateattributes(problem, {'Problem'}, {'nonempty'});

if ~ismember(class(problem), {'TerminalHitting'})
  error('Problem is not supported.');
end

end
