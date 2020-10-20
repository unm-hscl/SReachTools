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

surf(s, s, reshape(results.Pr(N-3, :), 100, 100), 'EdgeColor', 'none');

% width = 80;
% height = 115;
% gap = 10;
% 
% figure('Units', 'points', ...
%        'Position', [0, 0, 504, 150])
% 
% ax1 = subplot(1, 5, [1, 2], 'Units', 'points');
% data = reshape(Pr(N-3, :), 100, 100);
% ph = surf(ax1, s, s, data);
% ph.EdgeColor = 'none';
% caxis([0 1]);
% view([-30 45]);
% 
% % ax1.XGrid = 'off';
% % ax1.YGrid = 'off';
% % ax1.ZGrid = 'off';
% ax1.Color = 'none';
% 
% colorbar(ax1, 'off');
% ax1.Position = [40, 25, width*1.5, height];
% ax1.XLabel.Interpreter = 'latex';
% ax1.XLabel.String = '$\theta$';
% ax1.YLabel.Interpreter = 'latex';
% ax1.YLabel.String = '$\dot{\theta}$';
% ax1.Title.String = '(a)';
% set(ax1, 'FontSize', 8);
% 
% % ax2 = subplot(1, 5, 2, 'Units', 'points');
% % data = reshape(Pr(N-1, :), 100, 100);
% % ph = surf(ax2, s, s, data);
% % ph.EdgeColor = 'none';
% % caxis([0 1]);
% % view([0 90]);
% % 
% % colorbar(ax2, 'off');
% % % ax2.YAxis.Visible = 'off';
% % ax2.Position = [120, 25, width, height];
% % ax2.YLabel.Interpreter = 'latex';
% % ax2.YLabel.String = '$x_{2}$';
% % ax2.Title.String = '(b)';
% % set(ax2, 'FontSize', 8);
% 
% ax3 = subplot(1, 5, 3, 'Units', 'points');
% % data = reshape(Pr(1, :), 100, 100);
% data = reshape(Pr(N-3, :), 100, 100);
% ph = surf(ax3, s, s, data);
% ph.EdgeColor = 'none';
% caxis([0 1]);
% view([0 90]);
% 
% colorbar(ax3, 'off');
% % ax3.YAxis.Visible = 'off';
% ax3.Position = [120 + 1*(width + gap), 25, width, height];
% ax3.YLabel.Interpreter = 'latex';
% ax3.YLabel.String = '$\dot{\theta}$';
% ax3.Title.String = '(b)';
% 
% set(ax3, 'FontSize', 8);
% 
% %
% 
% ax4 = subplot(1, 5, 4, 'Units', 'points');
% data = reshape(Pr(N-3, :), 100, 100) + reshape(UB(N-3, :), 100, 100);
% ph = surf(ax4, s, s, data);
% ph.EdgeColor = 'none';
% caxis([0 1]);
% view([0 90]);
% 
% colorbar(ax4, 'off');
% ax4.YAxis.Visible = 'off';
% ax4.Position = [120 + 2*(width + gap), 25, width, height];
% ax4.YLabel.Interpreter = 'latex';
% ax4.YLabel.String = '$x_{2}$';
% ax4.Title.String = '(c)';
% 
% ax4.XLabel.Interpreter = 'latex';
% ax4.XLabel.String = '$\theta$';
% 
% set(ax4, 'FontSize', 8);
% 
% ax5 = subplot(1, 5, 5, 'Units', 'points');
% data = reshape(Pr(N-3, :), 100, 100) + reshape(LB(N-3, :), 100, 100);
% ph = surf(ax5, s, s, data);
% ph.EdgeColor = 'none';
% caxis([0 1]);
% view([0 90]);
% 
% colorbar(ax5);
% ax5.YAxis.Visible = 'off';
% ax5.Position = [120 + 3*(width + gap), 25, width, height];
% ax5.YLabel.Interpreter = 'latex';
% ax5.YLabel.String = '$x_{2}$';
% ax5.Title.String = '(d)';
% set(ax5, 'FontSize', 8);
