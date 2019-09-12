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
