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
            addParameter(p, 'ConstraintTube', srt.Tube.empty, ...
                @(x) isa(x, 'srt.Tube'));
            addParameter(p, 'TargetTube', srt.Tube.empty, ...
                @(x) isa(x, 'srt.Tube'));
            parse(p, varargin{:});

            obj = obj@srt.problems.Problem();

            obj.TargetTube = p.Results.TargetTube;
            obj.ConstraintTube = p.Results.ConstraintTube;

            % Ensure T in K?

            assert(length(obj.constraint_tube_) == ...
                    length(obj.target_tube_));

        end

        function tube = get.ConstraintTube(obj)
            tube = obj.constraint_tube_;
        end

        function tube = get.TargetTube(obj)
            tube = obj.target_tube_;
        end
    end

end
