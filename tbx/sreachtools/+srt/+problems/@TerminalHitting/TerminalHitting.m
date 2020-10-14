classdef TerminalHitting < srt.problems.Problem
% TERMINALHITTING Specifies a terminal-hitting time problem.
%
%   problem = TERMINALHITTING(N, K, T) creates a terminal-hitting time problem
%   object. The terminal-hitting time problem is defined as the probability of
%   reaching T at time k = N while remaining within K for all k < N.
%
%   See also: FirstHitting, MaximalHitting, Viability

    methods

        function obj = TerminalHitting(options)
            % TERMINALHITTING Construct an instance of the problem.

            arguments
                options.ConstraintTube (1, 1) srt.Tube
                options.TargetTube (1, 1) srt.Tube
            end

            obj = obj@srt.problems.Problem();

            % Ensure T in K?

            assert(length(options.ConstraintTube) == ...
                   length(options.TargetTube));

            obj.TargetTube = options.TargetTube;
            obj.ConstraintTube = options.ConstraintTube;

        end

    end

end
