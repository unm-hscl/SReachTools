%% Kernel Embeddings Example (Terminal-Hitting Time Problem)
% Kernel embeddings example showing the Terminal-hitting time problem
% for a Nonlinear Cart-Pole System
%
%%
% Specify the time horizon, the safe set $\mathcal{K}$, and the target set
% $\mathcal{T}$.
N = 1000;

K = srt.Tube(N, Polyhedron('lb', [-500 -500 -500 -500], 'ub', [500 500 500 500]));
T = srt.Tube(N, Polyhedron('lb', [-500 -500 -0.05 -500], 'ub', [500 500 0.05 500]));

problem = srt.problems.TerminalHitting('ConstraintTube', K, 'TargetTube', T);

%% System Definition
% Generate input/output samples for a nonlinear cart-pole system
%
% Load cart-pole data. 
load('..\data\CartPoleSamples_Nonlinear.mat');

% Select 5000 samples from the data.
indices = randperm(size(X, 2), 6000);
X = X(:, indices);
Y = Y(:, indices);
U = zeros(1, size(X, 2));

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
Xt = [
    zeros(1, size(Xt, 2));
    zeros(1, size(Xt, 2)); 
    Xt(1, :);
    Xt(2, :)
    ];
Ut = zeros(1, size(Xt, 2));

results = SReachPoint(problem, alg, sys, Xt, Ut);

%%
% View the results.

surf(s, s, reshape(results.Pr(N-1, :), 100, 100), 'EdgeColor', 'none');
