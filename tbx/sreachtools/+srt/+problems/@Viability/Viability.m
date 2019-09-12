classdef Viability < srt.problems.Problem
% VIABILITY Specifies a viability problem.
%
%   problem = VIABILITY(N, K) creates a viability problem object. The viability
%   problem is defined as the probability of remaining within K for all k < N.
%
%   See also: FirstHitting, MaximalHitting, TerminalHitting

    methods
        function obj = Viability(varargin)
            % VIABILITY Construct an instance of the problem.
            if nargin

                p = inputParser;
                addParameter(p, 'ConstraintTube', srt.Tube.empty);
                parse(p, varargin{:});

                obj.constraint_tube_ = p.Results.ConstraintTube;

            else
                error('Invalid problem definition.');
            end
        end
    end

end
