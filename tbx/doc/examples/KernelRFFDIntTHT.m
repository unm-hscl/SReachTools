%% Kernel Embeddings Example (Terminal-Hitting Time Problem)
% Kernel embeddings example showing the terminal-hitting time problem
% for a double integrator system.
%
%% Problem Definition
% Specify the time horizon, the safe set $\mathcal{K}$, and the target set
% $\mathcal{T}$.

N = 3;
K = srt.Tube(N, Polyhedron('lb', [-1; -1], 'ub', [1; 1]));
T = srt.Tube(N, Polyhedron('lb', [-0.5; -0.5], 'ub', [0.5; 0.5]));

prb = srt.problems.TerminalHitting('ConstraintTube', K, 'TargetTube', T);

%% System Definition
% Generate input/output samples for a double integrator system.
%
% $$x_{k+1} = A x_{k} + w_{k}, \quad w_{k} \sim \mathcal{N}(0, 0.01 I)$$
%

s = linspace(-1.1, 1.1, 50);
X = sampleunif(s, s);
U = zeros(1, size(X, 2));
W = 0.01.*randn(size(X));

A = [1, 0.25; 0, 1];
B = [0.03125; 0.25];

Y = A*X + B*U + W;

%%
% Create a sample-based stochastic system.

sys = srt.systems.SampledSystem('X', X, 'U', U, 'Y', Y);

%% Algorithm
% Initialize the algorithm.

alg = srt.algorithms.KernelEmbeddingsRFF('sigma', 1, 'lambda', 1, 'D', 500);

%%
% Call the algorithm.

s = linspace(-1, 1, 100);
X = sampleunif(s, s);
U = zeros(1, size(X, 2));

results = SReachPoint(prb, alg, sys, X, U);

%%
% View the results.

surf(s, s, reshape(results.Pr(1, :), 100, 100), 'EdgeColor', 'none');
