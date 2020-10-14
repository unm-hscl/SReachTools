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
%   High-level desc.   : Use Boole's inequality, Gaussian random vector,
%                        hyperbolic constraints-to-second order cone constraint
%                        reformulation, and piecewise linear approximation of
%                        the inverse of the standard normal cumulative density
%                        function to create a second-order cone program-based
%                        difference-of-convex optimization problem
%   Controller type    : A history-dependent affine controller that satisfies
%                        softened input constraints (controller satisfies the
%                        hard input bounds upto a user-specified probabilistic
%                        threshold)
%   Optimality         : Suboptimal affine controller for the
%                        underapproximation problem due to non-convexity
%                        established by the difference of convex formulation
%   Approximation      : Guaranteed underapproximation
%
%   This algorithm is part of the Stochastic Reachability Toolbox.
%   License for the use of this algorithm is given in
%   https://sreachtools.github.io/license/

    properties (Access = private)

        % Probabilistic relaxation of the hard input constraints
        max_input_viol_prob (1, 1) double {mustBePositive, ...
            mustBeLessThan(max_input_viol_prob, 1)} = 1E-2

        % Accuracy of piecewise-affine approximation of norminvcdf
        pwa_accuracy (1, 1) double {mustBePositive} = 1E-3

        % Difference-of-convex: Initialization of the slack multiplier
        tau_initial (1, 1) double {mustBePositive} = 1

        % Difference-of-convex: Scaling factor to the slack multiplier
        scaling_tau (1, 1) double {mustBePositive} = 2

        % Difference-of-convex: Max scaling factor
        tau_max (1, 1) double {mustBePositive} = 1E5

        % Difference-of-convex: Max iterations
        iter_max (1, 1) double {mustBePositive} = 200

        % Difference-of-convex: Exit condition tolerance for dc iterations
        dc_conv_tol (1, 1) double {mustBePositive} = 1E-4

        % Difference-of-convex: Slack tolerance requirements
        slack_tol (1, 1) double {mustBePositive} = 1E-8

    end

    methods
        function obj = ChanceAffine(varargin)
            % CHANCEAFFINE Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            p = inputParser;
            p.KeepUnmatched = true;

            addParameter(p, 'max_input_viol_prob', 1E-2);
            addParameter(p, 'pwa_accuracy', 1E-3);
            addParameter(p, 'tau_initial', 1);
            addParameter(p, 'scaling_tau', 2);
            addParameter(p, 'tau_max', 1E5);
            addParameter(p, 'iter_max', 200);
            addParameter(p, 'dc_conv_tol', 1E-4);
            addParameter(p, 'slack_tol', 1E-8);

            parse(p, varargin{:});

            obj.max_input_viol_prob = p.Results.max_input_viol_prob;
            obj.pwa_accuracy = p.Results.pwa_accuracy;
            obj.tau_initial = p.Results.tau_initial;
            obj.scaling_tau = p.Results.scaling_tau;
            obj.tau_max = p.Results.tau_max;
            obj.iter_max = p.Results.iter_max;
            obj.dc_conv_tol = p.Results.dc_conv_tol;
            obj.slack_tol = p.Results.slack_tol;


        end
    end

end
