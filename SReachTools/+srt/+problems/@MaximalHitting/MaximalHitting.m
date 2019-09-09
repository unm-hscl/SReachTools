classdef MaximalHitting < Problem
% MAXIMALHITTING Specifies a maximal-hitting time problem.
%
%   problem = MAXIMALHITTING(N, K, T) creates a maximal-hitting time problem
%   object. The maximal-hitting time problem is defined as...
%
%   See also: FirstHitting, TerminalHitting, Viability

  methods
    function obj = MaximalHitting(N, K, T, varargin)
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
        tf = ~obj.in_safe_set(varargin{:});
      end
    end
  end

end
