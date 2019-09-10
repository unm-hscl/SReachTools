function varargout = SReachPoint(prb, alg, sys, x0, varargin)
% SREACHPOINT Stochastic reachability.
%
%   SREACHPOINT(...) stochastic reachability for a point.
%
%   Example:
%       prob = srk.problems.FirstHittingTime();
%       alg = srk.algorithms.DynamicProgramming();
%       sys = srk.systems.StochasticLTISystem();
%       SReachPoint(prob, alg, sys);
%
%   See also: SReachForward, SReachSet
%
%   Copyright 2019 Adam Thorpe

p = inputParser;

valprob = @(arg) validateattributes(arg, {'Problem'}, {'nonempty'});
valalg  = @(arg) validateattributes(arg, {'Algorithm'}, {'nonempty'});
valsys  = @(arg) validateattributes(arg, {'StochasticSystem'}, {'nonempty'});
valvec  = @(arg) validateattributes(arg, {'numeric'}, {'column'});

addRequired(p, 'prb', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys', valsys);
addRequired(p, 'x0', valvec);
parse(p, prb, alg, sys, x0);

alg.compute_point(prb, sys, x0, varargin{:});

end
