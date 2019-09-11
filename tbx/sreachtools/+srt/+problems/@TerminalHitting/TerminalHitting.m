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
                addParameter(p, 'SafeSet', srt.Tube.empty);
                addParameter(p, 'TargetSet', srt.Tube.empty);
                parse(p, varargin{:});

                % Ensure T in K?

                obj.safe_tube_ = p.Results.SafeSet;
                obj.target_tube_ = p.Results.TargetSet;

                assert(length(obj.safe_tube_) == length(obj.target_tube_));

                obj.time_horizon_ = length(p.Results.SafeSet);
            else
%                 help(dbstack(1).name)
                error('Invalid problem definition.');
            end
        end
    end

    methods
        function tf = contains(obj, k, varargin)
            if k == obj.time_horizon_
                tube = obj.target_tube_(k);
                tf = tube.contains(varargin{:});
            else
                tube = obj.safe_tube_(k);
                tf = tube.contains(varargin{:});
            end
        end
    end

end
