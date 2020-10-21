classdef TerminalHitting < srt.problems.Problem
% TERMINALHITTING Specifies a terminal-hitting time problem.
%
%   problem = TERMINALHITTING(N, K, T) creates a terminal-hitting time problem
%   object. The terminal-hitting time problem is defined as the probability of
%   reaching T at time k = N while remaining within K for all k < N.
%
%   See also: FirstHitting, MaximalHitting, Viability

    methods

        function obj = TerminalHitting(varargin)
            % TERMINALHITTING Construct an instance of the problem.

            p = inputParser;
            p.KeepUnmatched = true;

            addParameter(p, 'ConstraintTube', srt.Tube.empty);
            addParameter(p, 'TargetTube', srt.Tube.empty);

            parse(p, varargin{:});

            obj = obj@srt.problems.Problem();

            % Ensure T in K?

            assert(length(p.Results.ConstraintTube) == ...
                   length(p.Results.TargetTube));

            obj.TargetTube = p.Results.TargetTube;
            obj.ConstraintTube = p.Results.ConstraintTube;

        end

    end

end
