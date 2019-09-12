function results = SReachPoint(prb, alg, sys, x0, varargin)
% SREACHPOINT Stochastic reachability.
%
%   SREACHPOINT(...) stochastic reachability for a point.
%
%   Example:
%       prob = srt.problems.FirstHittingTime();
%       alg = srt.algorithms.DynamicProgramming();
%       sys = srt.systems.StochasticLTISystem();
%       SReachPoint(prob, alg, sys);
%
%   See also: SReachForward, SReachSet
%
%   Copyright 2019 Adam Thorpe

import srt.*

p = inputParser;

valprob = @(arg) validateattributes(arg, {'srt.problems.Problem'}, {'nonempty'});
valalg  = @(arg) validateattributes(arg, {'srt.algorithms.Algorithm'}, {'nonempty'});
valsys  = @(arg) validateattributes(arg, {'srt.systems.StochasticSystem'}, {'nonempty'});
valvec  = @(arg) validateattributes(arg, {'numeric'}, {'nonempty'});

addRequired(p, 'prb', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys', valsys);
addRequired(p, 'x0', valvec);
parse(p, prb, alg, sys, x0);

results = alg.compute_point(prb, sys, x0, varargin{:});

end
