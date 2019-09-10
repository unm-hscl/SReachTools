function varargout = SReachSet(prb, alg, sys, varargin)
% SREACHSET Stochastic reachability.
%
%   SREACHSET(...) stochastic reachability for a set.
%
%   Example:
%       prob = srk.problems.FirstHittingTime();
%       alg = srk.algorithms.DynamicProgramming();
%       sys = srk.systems.StochasticLTISystem();
%       SReachSet(prob, alg, sys);
%
%   See also: SReachForward, SReachPoint
%
%   Copyright 2019 Adam Thorpe

p = inputParser;

valprob = @(arg) validateattributes(arg, {'Problem'}, {'nonempty'});
valalg  = @(arg) validateattributes(arg, {'Algorithm'}, {'nonempty'});
valsys  = @(arg) validateattributes(arg, {'StochasticSystem'}, {'nonempty'});

addRequired(p, 'prb', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys', valsys);
parse(p, prb, alg, sys);

alg.compute_set(prb, sys, varargin{:});

end
