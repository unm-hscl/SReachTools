function validateproblem(obj, problem)
% VALIDATEPROBLEM Validate problem.

supportedProblems = {'srt.problems.TerminalHitting'};

validateattributes(problem, supportedProblems, {'scalar', 'nonempty'});

end
