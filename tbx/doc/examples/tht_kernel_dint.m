%% Kernel Embeddings Example
% Kernel embeddings example 
%% 
% Specify the time horizon, the safe set $\mathcal{K}$, and the target set
% $\mathcal{T}$. 
N = 3;
K = srt.Tube(N, Polyhedron('lb', [-1.1; -1.1], 'ub', [1.1; 1.1]));
T = srt.Tube(N, Polyhedron('lb', [-1.1; -1.1], 'ub', [1.1; 1.1]));

% A = [1 0.25; 0 1];
% B = [0.03125; 0.25];
% U = Polyhedron('lb', -0.1, 'ub', 0.1);
%
% sys = StochasticLTISystem('A', A, 'B', B, ...
%                           'StateSpace', X, 'InputSpace', U);
%
% K = Polyhedron('lb', [-1; -1], 'ub', [1, 1]); % Safe set.
% T = Polyhedron('lb', [-1; -1], 'ub', [1, 1]); % Target set.
%
% problem = TerminalHitting(3, K, T);
%
% alg_dp = DynamicProgramming();
% alg_kr = KReach(SampleGenerator(sys));
%
% Reach('point', problem, alg_dp, sys, [0; 0]);
% Reach('point', problem, alg_kr, sys, [0; 0]);
