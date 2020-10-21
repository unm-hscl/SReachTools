function prob = getProbPolyhedron(obj, test_polyhedron, varargin)
% Compute the probability of a random vector lying in a polyhedron
% ====================================================================
% The probability computation is done via Monte-Carlo (quasi Monte-Carlo
% in Gaussian case), which is a guaranteed underapproximation (with a
% probabilistic guarantee). The guarantee in the general case comes
% via Hoeffding's inequality, and, in the Gaussian case, via the
% confidence interval.
%
% General case: A distribution-independent lower bound on the number
% of particles needed is obtained via Hoeffding's inequality.
%
% For \beta = exp(-2 * n_particles * \delta^2), Hoeffding's inequality
% states that
%
%               Prob{X - E[X] \geq \delta} \leq \beta,
%
% where, X is the empirical average of a collection of n_particles
% random variables bounded between [0,1], E[X] is the true mean, \delta
% is the desired_accuracy, \beta is the failure risk (the probability of
% the statement "X - E[X] \geq \delta" fails.
%
% In this function, we set the bounded random variables to be Bernoulli
% random variables that is 1 when a realization of the random vector
% is within the given test_polyhedron, and 0 otherwise.
% Consequently, E[X] simplifies to the probability of the
% realization to lie in the polyhedron.
%
% Given a failure risk \beta and desired accuracy \delta, we backcompute
% the number of particles as the following
%
%               n_particles = -ln(\beta) / (2 * delta^2)
%
% This enforces the condition that empirical probability estimate X does
% not overapproximate the true probability E[X] by more than delta.
%
% Gaussian case: We use Genz's algorithm to estimate the
% probability. We increase the number of particles in the powers of
% 10, till error estimate is within the desired tolerance.
%
% In both of the cases, to remain conservative in our estimate, we
% subtract the desired_accuracy from the emperical probability estimate.
% In other words, we use the lower bound of the confidence interval.
% Thus, the probability estimate is guaranteed to be an
% underapproximation.
%
% Inputs:
% -------
%   obj             - RandomVector object
%   test_polyhedron - Polyhedron object (polytope whose probability of
%                     occurrence is of interest)
%   desired_accuracy- [Optional] Maximum absolute deviation from the
%                     true probability estimate [Default: 1e-2]
%
% Outputs:
% --------
%   covar           - Probability of the random vector lying in the
%                     given polytope
% Notes:
% ------
% * Due to the inverse-square dependence on the desired_accuracy, we
%   impose a hard lower bound of 1e-2 on the desired_accuracy. This
%   leads to the requirement of 2e5 particles.
% * We set the default failure risk as 2e-15.
% * Ill-formed (mean/cov has Inf/NaN) Gaussian random vectors return
%   zero probability.
% * We compute the floor of the probability value based on the
%   desired_accuracy to obtain a preferable lower bound on the
%   probability.
%
% ====================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%

dim = size(obj.Mean(), 1);

failure_risk = 2e-15;
gauss_scale_factor = 10;
if ~isempty(varargin)
    desired_accuracy = varargin{1};
    validateattributes(desired_accuracy, {'numeric'},...
        {'scalar','>=', 1e-2, '<=', 1});
    if length(varargin) > 1
        throwAsCaller(SrtInvalidArgsError('Too many inputs'));
    end
else
    desired_accuracy = 1e-2;
end

% Check if we got a MPT's Polyhedron object of correct dimensions
validateattributes(test_polyhedron, {'Polyhedron'}, {'nonempty'},...
    'RandomVector/getProbPolyhedron','test_polyhedron');
if test_polyhedron.Dim ~= dim
    throwAsCaller(SrtInvalidArgsError(['Mismatch in polytope ',...
        'dimensions and random vector dimensions']));
end

if any(isinf(abs(obj.Mean()))) || ...
        any(isnan(abs(obj.Mean()))) || ...
        any(any(isinf(abs(obj.Sigma())))) || ...
        any(any(isnan(abs(obj.Sigma()))))
    % Return zero probability if ill-formed random vector
    temp_probability = 0;
else
    % Construct the half-space representation for qscmvnv
    qscmvnv_lb = repmat(-Inf, ...
        [size(test_polyhedron.A, 1),1]);
    qscmvnv_coeff_matrix = test_polyhedron.A;
    % We have the polytope given by Ax <= b, but qscmvnv
    % expects a zero mean Gaussian. Hence, define x->x+mean
    qscmvnv_ub = test_polyhedron.b - ...
        test_polyhedron.A*obj.Mean();
    % Set up for an iterative call of Genz's particles;
    % Increase particles till the error estimate is within
    % threshold
    n_particles = 10;
    est_3sig_CI = 1;
    while est_3sig_CI > desired_accuracy
        try
            % Call Genz's algorithm;
            [temp_probability, est_3sig_CI] = qscmvnv( ...
                n_particles,obj.Sigma(), qscmvnv_lb, ...
                qscmvnv_coeff_matrix, qscmvnv_ub);
        catch ME
            disp(ME);
            throw(SrtDevError(['Error in qscmvnv (', ...
                'Quadrature of multivariate Gaussian)']));
        end
        n_particles = n_particles * gauss_scale_factor;
    end
end

% Subtract the desired_accuracy to remain an underapproximation
temp_probability = temp_probability - desired_accuracy;

if temp_probability > 0
    % Rounding DOWN the integral to the desired accuracy
    prob = floor(temp_probability/desired_accuracy) * ...
        desired_accuracy;
else
    prob = 0;
end

end
