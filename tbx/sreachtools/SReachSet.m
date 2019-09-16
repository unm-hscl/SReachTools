function results = SReachSet(prb, alg, sys, varargin)
% SREACHSET Stochastic reachability.
%
%   SREACHSET(...) stochastic reachability for a set.
%
%   Example:
%       prob = srt.problems.FirstHittingTime();
%       alg = srt.algorithms.DynamicProgramming();
%       sys = srt.systems.StochasticLTISystem();
%       SReachSet(prob, alg, sys);
%
%   See also: SReachForward, SReachPoint
%
%   Copyright 2019 Adam Thorpe

import srt.*

p = inputParser;

valprob = @(arg) isa(arg, 'srt.problems.Problem');
valalg  = @(arg) isa(arg, 'srt.algorithms.Algorithm');
valsys  = @(arg) isa(arg, 'srt.systems.StochasticSystem');

addRequired(p, 'prb', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys', valsys);
parse(p, prb, alg, sys);

results = alg.compute_set(prb, sys, varargin{:});

end
