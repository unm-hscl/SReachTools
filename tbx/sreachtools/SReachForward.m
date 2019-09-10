function varargout = SReachForward(problem, alg, sys, varargin)
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

p = inputParser;

valprob = @(arg) validateattributes(arg, {'Problem'}, {'nonempty'});
valalg  = @(arg) validateattributes(arg, {'Algorithm'}, {'nonempty'});

addRequired(p, 'problem', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys', validatesystem(sys));
parse(p, problem, alg, sys);

alg.validatearguments(problem, sys);
alg.compute_fwd(problem, sys, varargin{:});

end
