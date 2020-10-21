classdef FirstHitting < srt.problems.Problem
% FIRSTHITTING Specifies a first-hitting time problem.
%
%   problem = FIRSTHITTING(N, K, T) creates a first-hitting time problem object.
%   The first-hitting time problem is defined as the probability of reaching T
%   at time j <= N while remaining within K for all k < j.
%
%   See also: MaximalHitting, TerminalHitting, Viability

    methods

        function obj = FirstHitting(varargin)
            % FIRSTHITTING Construct an instance of the problem.

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
