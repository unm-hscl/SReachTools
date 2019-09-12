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
            if nargin

                p = inputParser;
                addParameter(p, 'ConstraintTube', srt.Tube.empty);
                addParameter(p, 'TargetTube', srt.Tube.empty);
                parse(p, varargin{:});

                % Ensure T in K?

                obj.constraint_tube_ = p.Results.ConstraintTube;
                obj.target_tube_ = p.Results.TargetTube;

                assert(length(obj.constraint_tube_) == ...
                       length(obj.target_tube_));

            else
                error('Invalid problem definition.');
            end
        end
    end

end
