classdef ChanceAffineUniform < srt.algorithms.Algorithm
% CHANCEAFFINEUNIFORM Chance-constrained affine uniform.
%
%   Solve the problem of stochastic reachability of a target tube (a lower bound
%   on the maximal reach probability and an affine controller synthesis) using
%   chance-constrained optimization and uniform risk allocation
%
%   SReachPointCcAuniform implements the chance-constrained underapproximation
%   to the problem of stochastic reachability of a target tube to construct an
%   affine controller. This technique is inspired from Algorithms 1 and 2 of
%
%   M. Vitus and C. Tomlin, "On feedback design and risk allocation in chance
%   constrained control", In Proc. Conf. Dec. & Ctrl., 2011.
%
%   In contrast to their original algorithm, we have a chance constraint on the
%   input and the state. Further, the lower bound on the reachability (state
%   constraint) depends on how high the input chance constraint satisfaction
%   probability is. Therefore, we perform two levels of bisection --- one to
%   maximize the probability of constraint satisfaction for the state, and the
%   other to meet the chance constraint on the input. However, to save time, we
%   check only for feasibility in the input bisection.
%
%   Subsequently, the obtained solution is discounted for input constraint
%   violation using Theorem 1 of
%
%   A. Vinod and M. Oishi. Affine controller synthesis for stochastic
%   reachability via difference of convex programming. In Proc. Conf. Dec. &
%   Ctrl., 2019. (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
%
%   This algorithm is part of the Stochastic Reachability Toolbox.
%   License for the use of this algorithm is given in
%   https://sreachtools.github.io/license/

    methods
        function obj = ChanceAffineUniform(varargin)
            % CHANCEAFFINEUNIFORM Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

        end
    end

end
