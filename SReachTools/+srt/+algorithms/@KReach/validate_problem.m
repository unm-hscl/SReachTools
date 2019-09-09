function validate_problem(obj, problem)
% VALIDATE_PROBLEM Checks if the problem is valid for the algorithm.

validateattributes(problem, {'Problem'}, {'nonempty'});

if ~ismember(class(problem), {'TerminalHitting', 'FirstHitting'})
  error('Problem is not supported.');
end

end
