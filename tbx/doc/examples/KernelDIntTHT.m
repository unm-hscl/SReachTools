%% Kernel Embeddings Example (Terminal-Hitting Time Problem)
% Kernel embeddings example showing the terminal-hitting time problem
% for a double integrator system.
%
%% Problem Definition
% Specify the time horizon, the safe set $\mathcal{K}$, and the target set
% $\mathcal{T}$.

N = 1000;
% K = srt.Tube(N, Polyhedron('lb', [-1; -1], 'ub', [1; 1]));
% T = srt.Tube(N, Polyhedron('lb', [-0.5; -0.5], 'ub', [0.5; 0.5]));
K = srt.Tube(N, @(x) -0.05 <= x(2) & x(2) <= 0.05);
T = srt.Tube(N, @(x) -pi/6 <= x(2) & x(2) <= pi/6);

problem = srt.problems.TerminalHitting('ConstraintTube', K, 'TargetTube', T);

%% System Definition
% Generate input/output samples for a double integrator system.
%
% $$x_{k+1} = A x_{k} + w_{k}, \quad w_{k} \sim \mathcal{N}(0, 0.01 I)$$
%

% s = linspace(-1.1, 1.1, 5);
% X = sampleunif(s, s);
% U = zeros(1, size(X, 2));
% W = 0.01.*randn(size(X));
%
% A = [1, 0.25; 0, 1];
% B = [0.03125; 0.25];
%
% Y = A*X + B*U + W;

X = [data1 data2 data3 data4];

%%
% Create a sample-based stochastic system.

sys = srt.systems.SampledSystem('X', X, 'U', U, 'Y', Y);

%% Algorithm
% Initialize the algorithm.

alg = srt.algorithms.KernelEmbeddings('sigma', 0.1, 'lambda', 1);

%%
% Call the algorithm.

s = linspace(-1, 1, 100);
Xt = sampleunif(s, s);
Ut = zeros(1, size(Xt, 2));

results = SReachPoint(problem, alg, sys, Xt, Ut);

%%
% View the results.

surf(s, s, reshape(results.Pr(1, :), 100, 100), 'EdgeColor', 'none');

% xlabel('$$\theta$$', 'Interpreter', 'latex')
% ylabel('$$\dot{\theta}$$', 'Interpreter', 'latex')
% zlabel('$$r_{x_{0}}^{\pi}(\mathcal{K}, \mathcal{T})$$', 'Interpreter', 'latex')
