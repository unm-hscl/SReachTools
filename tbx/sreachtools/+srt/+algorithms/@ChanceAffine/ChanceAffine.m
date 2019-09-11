classdef ChanceAffine < srt.algorithms.Algorithm
% CHANCEAFFINE Chance-constrained affine.
%
%   Solve the problem of stochastic reachability of a target tube (a lower bound
%   on the maximal reach probability and an affine controller synthesis) using
%   chance-constrained optimization and difference of convex programming.
%
%   ChanceAffine implements the chance-constrained underapproximation to the
%   problem of stochastic reachability of a target tube to construct an affine
%   controller. This technique is discussed in detail in the paper,
%
%   A. Vinod and M. Oishi. Affine controller synthesis for stochastic
%   reachability via difference of convex programming. In Proc. Conf. Dec. &
%   Ctrl., 2019. (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
%
%    High-level desc.   : Use Boole's inequality, Gaussian random vector,
%                         hyperbolic constraints-to-second order cone constraint
%                         reformulation, and piecewise linear approximation of
%                         the inverse of the standard normal cumulative density
%                         function to create a second-order cone program-based
%                         difference-of-convex optimization problem
%    Controller type    : A history-dependent affine controller that satisfies
%                         softened input constraints (controller satisfies the
%                         hard input bounds upto a user-specified probabilistic
%                         threshold)
%    Optimality         : Suboptimal affine controller for the
%                         underapproximation problem due to non-convexity
%                         established by the difference of convex formulation
%    Approximation      : Guaranteed underapproximation
%
%   This algorithm is part of the Stochastic Reachability Toolbox.
%   License for the use of this algorithm is given in
%   https://sreachtools.github.io/license/

    methods
        function obj = ChanceAffine(varargin)
            % CHANCEAFFINE Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

        end
    end

end
