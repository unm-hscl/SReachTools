%% Controller synthesis using |SReachPoint| for a Dubin's vehicle
%   Chance-constrained approach that uses risk allocation and piecewise-affine
%   approximations to formulate a difference-of-convex program to synthesize a
%   closed-loop (affine disturbance feedback) controller. The controller
%   synthesis is done by solving a series of second-order cone programs. (See
%   <http://hscl.unm.edu/affinecontrollersynthesis Vinod and Oishi, Conference
%   in Decision and Control, 2019 (submitted)>)

%% Dubin's vehicle dynamics
% We consider a Dubin's vehicle with known turning rate sequence
% $\overline{\omega} = {[\omega_0\ \omega_1\ \ldots\ \omega_{T-1}]}^\top
% \in R^T$, with additive Gaussian disturbance. Specifically, we consider
% $T=25$ and set the turning rate to $\omega_\mathrm{max}=\frac{\pi}{T_s T}$ for
% the first half of the time interval, and $-\omega_\mathrm{max}$ for the rest
% of the time interval. Here, $T_s=0.1$ is the sampling time. The resulting
% dynamics are,
%
% $$x_{k+1} = x_k + T_s \cos\left(\theta_0 + \sum_{i=1}^{k-1}
% \omega_i T_s\right) v_k + \eta^x_k$$
%
% $$y_{k+1} = y_k + T_s \sin\left(\theta_0 + \sum_{i=1}^{k-1}
% \omega_i T_s\right) v_k + \eta^y_k$$
%
% where $x,y$ are the positions (state) of the Dubin's vehicle in $\mathrm{x}$-
% and $\mathrm{y}$- axes, $v_k$ is the velocity of the vehicle (input),
% $\eta^{(\cdot)}_k$ is the additive Gaussian disturbance affecting the
% dynamics, and $\theta_0$ is the initial heading direction. We define the
% disturbance as ${[\eta^x_k\ \eta^y_k]}^\top\sim \mathcal{N}({[0\ 0]}^\top,
% 0.005 I_2)$.

n_mcarlo_sims = 1e5;                        % Monte-Carlo simulation particles
n_mcarlo_sims_affine = 1e3;                 % For affine controllers

N = 50;                          % Time Horizon
Ts = 0.1;                        % Sampling time

umax = 10;
InputSpace = Polyhedron('lb', 0, 'ub', umax);

prob_thresh = 0.9;

sys = srtDubinsCarModel(N, Ts, ...
    'DisturbanceType', 'affine', ...
    'InputSpace', InputSpace);

%% Target tube definition
% We define the target tube to be a collection of time-varying boxes
% $\{\mathcal{T}_k\}_{k=0}^N$ where $N$ is the time horizon.
%
% In this problem, we define $\mathcal{T}_k$ to be centered about the nominal
% trajectory with fixed velocity of $u_\mathrm{max} * 3/2$ (faster than the
% maximum velocity allowed) and the heading angle sequence with $\pi/2$ removed.
% The half-length of these boxes decay exponentially with a time constant which
% is $N/2$.

v_nominal = umax * 2/3;                 % Nominal trajectory's heading velocity
% Construct the nominal trajectory
[~,H,~] = sys.getConcatMats(N);
center_box_X = [zeros(2,1);H * (v_nominal * ones(N,1))];
center_box = reshape(center_box_X,2,[]);
% Box sizes
box_halflength_at_0 = 4;                % Box half-length at t=0
time_const = 1/2*N;          % Time constant characterize the
                                        % exponentially decaying box half-length


% Target tube definition as well as plotting
target_tube_cell = cell(N + 1,1); % Vector to store target sets
figure(100);clf;hold on

for itt = 0:N
    % Define the target set at time itt
    target_tube_cell{itt+1} = Polyhedron(...
        'lb',center_box(:, itt+1) -box_halflength_at_0*exp(- itt/time_const),...
        'ub', center_box(:, itt+1) + box_halflength_at_0*exp(- itt/time_const));
end

axis equal

% Target tube definition
target_tube = srt.Tube(target_tube_cell{:});

h_target_tube = plot(target_tube);

h_nominal_traj = scatter(center_box(1,:), center_box(2,:), 50,'ks','filled');
h_vec = [h_target_tube, h_nominal_traj];
legend_cell = {'Target tube', 'Nominal trajectory'};
legend(h_vec, legend_cell, 'Location','EastOutside', 'interpreter','latex');
xlabel('x');
ylabel('y');
axis equal
box on;
grid on;
drawnow;

% Initial state.
init_state = [-1.5;1.5];

%% Quantities needed to compute the optimal mean trajectory
% We first compute the dynamics of the concatenated state vector $X = Z x_0
% + H U + G W$, and compute the concatentated random vector $W$ and its mean.
[Z,H,G] = sys.getConcatMats(N);
% Compute the mean trajectory of the concatenated disturbance vector
muW_gauss = sys.Disturbance.concat(N).Mean();

%% |SReachPoint|: |chance-open|
% This method is discussed in <http://hscl.unm.edu/affinecontrollersynthesis
% Vinod and Oishi, Hybrid Systems: Computation and Control, 2019 (submitted)>.
% It was introduced for stochastic reachability in
% <http://doi.org/10.1109/CDC.2013.6760626 Lesser et. al., Conference on
% Decision and Control, 2013>.
%
% This approach implements the chance-constrained approach to compute an optimal
% open-loop controller. It uses risk allocation and piecewise-affine
% overapproximation of the inverse normal cumulative density function to
% formulate a linear program for this purpose. Naturally, this is one of the
% fastest ways to compute an open-loop controller and an underapproximative
% probabilistic guarantee of safety. However, due to the use of Boole's
% inequality for risk allocation, it provides a conservative estimate of safety
% using the open-loop controller.

fprintf('\n\nSReachPoint with chance-open\n');

prob = srt.problems.TerminalHitting( ...
    'TargetTube', target_tube, ...
    'ConstraintTube', target_tube);

alg = srt.algorithms.ChanceOpen('pwa_accuracy', 1E-3);

timerVal = tic;

results = SReachPoint(prob, alg, sys, init_state);

elapsed_time_particle_gauss = toc(timerVal);

if prob_chance_open_gauss
    % Optimal mean trajectory construction
    % mean_X = Z * x_0 + H * U + G * \mu_W
    opt_mean_X_chance_open_gauss = Z * init_state_chance_open_gauss + ...
        H * opt_input_vec_chance_open_gauss + G * muW_gauss;
    opt_mean_traj_chance_open_gauss = reshape( ...
        opt_mean_X_chance_open_gauss, sys_gauss.state_dim,[]);
    % Check via Monte-Carlo simulation
    concat_state_realization = generateMonteCarloSims(n_mcarlo_sims, ...
        sys_gauss, init_state_chance_open_gauss, time_horizon,...
        opt_input_vec_chance_open_gauss);
    mcarlo_result = target_tube.contains(concat_state_realization);
    simulated_prob_chance_open_gauss = sum(mcarlo_result)/n_mcarlo_sims;
else
    simulated_prob_chance_open_gauss = NaN;
end

fprintf('SReachPoint underapprox. prob: %1.2f | Simulated prob: %1.2f\n',...
    prob_chance_open_gauss, simulated_prob_chance_open_gauss);
fprintf('Computation time: %1.3f\n', elapsed_time_chance_open_gauss);
