function validateproblem(obj, problem)
% VALIDATEPROBLEM Checks if the problem is valid for the algorithm.

validateattributes(problem, {'srt.problems.Problem'}, {'nonempty'});

supportedProblems = {'srt.problems.TerminalHitting', ...
                     'srt.problems.FirstHitting'};

if ~ismember(class(problem), supportedProblems)
  error('Problem is not supported.');
end

end
