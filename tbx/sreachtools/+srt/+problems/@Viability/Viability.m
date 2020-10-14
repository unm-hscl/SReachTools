classdef Viability < srt.problems.Problem
% VIABILITY Specifies a viability problem.
%
%   problem = VIABILITY(N, K) creates a viability problem object. The viability
%   problem is defined as the probability of remaining within K for all k < N.
%
%   See also: FirstHitting, MaximalHitting, TerminalHitting

    methods

        function obj = Viability(options)
            % VIABILITY Construct an instance of the problem.

            arguments
                options.ConstraintTube (1, 1) srt.Tube
            end

            obj = obj@srt.problems.Problem();

            obj.ConstraintTube = options.ConstraintTube;

        end

    end

end
