classdef FirstHitting < Problem
% FIRSTHITTING Specifies a first-hitting time problem.
%
%   problem = FIRSTHITTING(N, K, T) creates a first-hitting time problem object.
%   The first-hitting time problem is defined as the probability of reaching T
%   at time j <= N while remaining within K for all k < j.
%
%   See also: MaximalHitting, TerminalHitting, Viability

    methods
        function obj = FirstHitting(N, K, T, varargin)
            p = inputParser;
            addRequired(p, 'N');
            addRequired(p, 'K');
            addRequired(p, 'T');
            parse(p, N, K, T, varargin{:});

            obj.time_horizon_ = N;
            obj.safe_set_ = K;
            obj.target_set_ = T;
        end
    end

    methods
        function tf = contains(obj, k, varargin)
            if k == obj.time_horizon_
                tf = ~obj.in_target_set(varargin{:});
            else
                tf = ~obj.in_safe_set(varargin{:}) | ...
                      obj.in_target_set(varargin{:});
            end
        end
    end

end
