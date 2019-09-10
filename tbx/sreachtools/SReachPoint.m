function varargout = SReachPoint(problem, alg, sys, varargin)
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

addRequired(p, 'problem', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys', validatesystem(sys));
parse(p, problem, alg, sys);

alg.validatearguments(problem, sys);
alg.compute_point(problem, sys, varargin{:});

end
