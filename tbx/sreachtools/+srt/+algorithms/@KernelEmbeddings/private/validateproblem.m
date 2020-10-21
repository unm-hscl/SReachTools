function validateproblem(obj, problem)
% VALIDATEPROBLEM Checks if the problem is valid for the algorithm.

supportedProblems = {'srt.problems.FirstHitting', ...
                     'srt.problems.TerminalHitting', ...
                     'srt.problems.Viability'};

validateattributes(problem, supportedProblems, {'scalar', 'nonempty'});

end
