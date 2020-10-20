%% Kernel Embeddings Example (Terminal-Hitting Time Problem)
% Kernel embeddings example showing the Terminal-hitting time problem
% for a Linear Cart-Pole System
%
%%
% Specify the time horizon, the safe set $\mathcal{K}$, and the target set
% $\mathcal{T}$.
N = 20;

K = srt.Tube(N, Polyhedron('lb', [-0.7 -1 -pi\6 -500], 'ub', [0.7 1 pi\6 500]));
T = srt.Tube(N, Polyhedron('lb', [-500 -500 -0.05 -500], 'ub', [500 500 0.05 500]));

problem = srt.problems.TerminalHitting('ConstraintTube', K, 'TargetTube', T);

%% System Definition
% Generate input/output samples for a linear cart-pole system
%
% Load cart-pole data. 
load('..\data\CartPoleSamples_Linearized.mat');

% Select 5000 samples from the data.
indices = randperm(size(X, 2), 5000);
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

width = 80;
height = 100;

figure('Units', 'points', ...
       'Position', [0, 0, 200, 150])

ax = axes; %subplot(1, 5, [1, 2], 'Units', 'points');
ax.Units = 'points';
data = reshape(results.Pr(N-3, :), 100, 100);
ph = surf(ax, s, s, data);
ph.EdgeColor = 'none';
caxis([0 1]);
view(2);
ax.Color = 'none';

colorbar(ax, 'off');
ax.Position = [40, 25, width*1.5, height];
ax.XLabel.Interpreter = 'latex';
ax.XLabel.String = '$\theta$';
ax.YLabel.Interpreter = 'latex';
ax.YLabel.String = '$\dot{\theta}$';
ax.Title.String = '(a)';
set(ax, 'FontSize', 8);

saveas(gcf,'.\images\KernelLinearCartPoleTHT','eps')
saveas(gcf,'.\images\KernelLinearCartPoleTHT','fig')
