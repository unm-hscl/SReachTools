classdef Viability < Problem
% VIABILITY Specifies a viability problem.
%
%   problem = VIABILITY(N, K) creates a viability problem object. The viability
%   problem is defined as the probability of remaining within K for all k < N.
%
%   See also: FirstHitting, MaximalHitting, TerminalHitting

  methods
    function obj = Viability(N, K, varargin)
      p = inputParser;
      addRequired(p, 'N');
      addRequired(p, 'K');
      parse(p, N, K, varargin{:});

      obj.time_horizon_ = N;
      obj.safe_set_ = K;
    end
  end

  methods
    function tf = contains(obj, k, varargin)
      tf = ~obj.in_safe_set(varargin{:});
    end
  end

end
