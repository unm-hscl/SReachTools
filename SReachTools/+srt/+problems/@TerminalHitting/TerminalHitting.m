classdef TerminalHitting < Problem
% TERMINALHITTING Specifies a terminal-hitting time problem.
%
%   problem = TERMINALHITTING(N, K, T) creates a terminal-hitting time problem
%   object. The terminal-hitting time problem is defined as the probability of
%   reaching T at time k = N while remaining within K for all k < N.
%
%   See also: FirstHitting, MaximalHitting, Viability

  methods
    function obj = TerminalHitting(N, K, T, varargin)
      p = inputParser;
      addRequired(p, 'N');
      addRequired(p, 'K');
      addRequired(p, 'T');
      parse(p, N, K, T, varargin{:});

      % Ensure T in K?

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
        tf = ~obj.in_safe_set(varargin{:});
      end
    end
  end

end
