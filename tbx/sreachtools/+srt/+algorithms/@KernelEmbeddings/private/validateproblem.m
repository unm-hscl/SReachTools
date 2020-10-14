function validateproblem(obj, problem)
% VALIDATEPROBLEM Checks if the problem is valid for the algorithm.

arguments
    obj
    problem (1, 1) srt.problems.Problem {mustBeNonempty}
end

supportedProblems = {'srt.problems.FirstHitting', ...
                     'srt.problems.TerminalHitting', ...
                     'srt.problems.Viability'};

if ~ismember(class(problem), supportedProblems)
  error('Problem is not supported.');
end

end
