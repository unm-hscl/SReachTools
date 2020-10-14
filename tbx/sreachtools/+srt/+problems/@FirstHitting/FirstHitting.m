classdef FirstHitting < srt.problems.Problem
% FIRSTHITTING Specifies a first-hitting time problem.
%
%   problem = FIRSTHITTING(N, K, T) creates a first-hitting time problem object.
%   The first-hitting time problem is defined as the probability of reaching T
%   at time j <= N while remaining within K for all k < j.
%
%   See also: MaximalHitting, TerminalHitting, Viability

    methods

        function obj = FirstHitting(options)
            % FIRSTHITTING Construct an instance of the problem.

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
