function results = SReachSet(prb, alg, sys, lev, varargin)
% SREACHSET Stochastic reachability.
%
%   SREACHSET(...) stochastic reachability for a set.
%
%   Example:
%       prob = srt.problems.FirstHittingTime();
%       alg = srt.algorithms.DynamicProgramming();
%       sys = srt.systems.StochasticLTISystem();
%       lev = 0.8
%       SReachSet(prob, alg, sys, lev);
%
%   See also: SReachForward, SReachPoint
%
%   Copyright 2019 Adam Thorpe

import srt.*

p = inputParser;

vallev   = @(arg) validateattributes(arg, {'numeric'}, ...
    {'scalar', 'nonnegative', '<=', 1});
valprob = @(arg) isa(arg, 'srt.problems.Problem');
valalg  = @(arg) isa(arg, 'srt.algorithms.Algorithm');
valsys  = @(arg) isa(arg, 'srt.systems.StochasticSystem');

addRequired(p, 'prb', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys', valsys);
addRequired(p, 'lev', vallev);
parse(p, prb, alg, sys, lev);

results = alg.compute_set(prb, sys, lev, varargin{:});

end
