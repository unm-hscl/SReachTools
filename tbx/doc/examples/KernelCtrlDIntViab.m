%% Kernel Embeddings Example (Viability Problem)
% Kernel embeddings example (controller synth) showing the viability problem
% for a double integrator system.
%
%%
% Specify the time horizon, the safe set $\mathcal{K}$, and the target set
% $\mathcal{T}$.

N = 3;
K = srt.Tube(N, Polyhedron('lb', [-1; -1], 'ub', [1; 1]));

problem = srt.problems.Viability('ConstraintTube', K);

%% System Definition
% Generate input/output samples for a double integrator system.
%
% $$x_{k+1} = A x_{k} + w_{k}, \quad w_{k} \sim \mathcal{N}(0, 0.01 I)$$
%

X = -1.1 + (2.2).*rand(2, 500);
U = -0.1 + (0.2).*rand(1, 500);

W = 0.01.*randn(size(X));

A = [1, 0.25; 0, 1];
B = [0.03125; 0.25];

Y = A*X + B*U + W;

%%
% Create a sample-based stochastic system.

sys = srt.systems.SampledSystem('X', X, 'U', U, 'Y', Y);

%% Create test points
%

s = linspace(-1, 1, 50);
X = sampleunif(s, s);
U = zeros(1, size(X, 2));

%% Compute safety probabilities
% (uncontrolled)

alg = srt.algorithms.KernelEmbeddings('sigma', 0.1, 'lambda', 1, ...
                                      'verbose', 1);

results_uncontrolled = SReachPoint(problem, alg, sys, X, U);

%% Compute safety probabilities
% (controlled)

alg = srt.algorithms.KernelEmbeddingsControl('sigma', 0.1, 'lambda', 1, ...
                                             'verbose', 2);

results_controlled = SReachPoint(problem, alg, sys, X);

%% Plot the Results
%

figure
surf(s, s, reshape(results_uncontrolled.Pr(1, :), 50, 50), ...
     'EdgeColor', 'none');
view([0 90]);

figure
surf(s, s, reshape(results_controlled.Pr(1, :), 50, 50), ...
     'EdgeColor', 'none');
view([0 90]);

d = results_controlled.Pr(1, :) - results_uncontrolled.Pr(1, :);

figure
surf(s, s, reshape(d, 50, 50), ...
     'EdgeColor', 'none');
caxis([0 1]);zlim([0 1]);
view([0 90]);
