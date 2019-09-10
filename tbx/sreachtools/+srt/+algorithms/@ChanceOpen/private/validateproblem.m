function validateproblem(obj, prob)
% VALIDATEPROBLEM Validate problem.

validateattributes(problem, {'Problem'}, {'nonempty'});

if ~ismember(class(problem), {'TerminalHitting'})
  error('Problem is not supported.');
end

% if ~(strcmpi(options.prob_str, 'term') &&...
%         strcmpi(options.method_str, 'chance-open'))
%     throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
% end

end
