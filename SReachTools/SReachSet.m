function varargout = SReachSet(problem, alg, sys, varargin)
% SREACHSET Stochastic reachability.
%
%   SREACHSET(...) stochastic reachability for a set.

p = inputParser;

valprob = @(arg) validateattributes(arg, {'Problem', 'char'}, {'nonempty'});
valalg  = @(arg) validateattributes(arg, {'Algorithm'}, {'nonempty'});

addRequired(p, 'problem', valprob);
addRequired(p, 'alg', valalg);
addRequired(p, 'sys');
parse(p, problem, alg, sys);

if isa(problem, 'char')
  switch problem
    case 'term'
      % Do something.
    case 'first'
      % Do something.
    case 'max'
      % Do something.
    otherwise
      error('Not a valid problem.');
  end
end

alg.validate_arguments(problem, sys);
alg.compute_set(problem, sys, varargin{:});

end
