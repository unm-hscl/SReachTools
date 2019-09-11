function results = SReachForward(prb, alg, sys, varargin)
% SREACHFORWARD Forward stochastic reachability.
%
%   SREACHFORWARD(...) forward stochastic reachability.
%
%   Example:
%       prob = srk.problems.FirstHittingTime();
%       alg = srk.algorithms.DynamicProgramming();
%       sys = srk.systems.StochasticLTISystem();
%       SReachForward(prob, alg, sys);
%
%   See also: SReachPoint, SReachSet
%
%   Copyright 2019 Adam Thorpe

import srt.*

p = inputParser;

valprob = @(arg) validateattributes(arg, {'Problem'}, {'nonempty'});
valalg  = @(arg) validateattributes(arg, {'Algorithm'}, {'nonempty'});
valsys  = @(arg) validateattributes(arg, {'StochasticSystem'}, {'nonempty'});

addRequired(p, 'prb', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys', valsys);
parse(p, prb, alg, sys);

results = alg.compute_fwd(prb, sys, varargin{:});

end
