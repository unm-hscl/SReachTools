function results = SReachForward(prb, alg, sys, varargin)
% SREACHFORWARD Forward stochastic reachability.
%
%   SREACHFORWARD(...) forward stochastic reachability.
%
%   Example:
%       prob = srt.problems.FirstHittingTime();
%       alg = srt.algorithms.DynamicProgramming();
%       sys = srt.systems.StochasticLTISystem();
%       SReachForward(prob, alg, sys);
%
%   See also: SReachPoint, SReachSet
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

results = alg.compute_fwd(prb, sys, varargin{:});

end
