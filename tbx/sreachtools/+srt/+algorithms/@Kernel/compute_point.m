function result = compute_point(obj, problem, sys, x0, varargin)

obj.validate_dependencies();
p = inputParser;
addRequired(p, 'problem', @obj.validate_problem);
addRequired(p, 'sys', @obj.validate_system);
addRequired(p, 'x0');
parse(p, problem, sys, x0, varargin{:});

% Constants
N = problem.TimeHorizon;

mt = size(x0, 2);

M = obj.num_samples_;

% Compute weight matrix.
Gx = obj.compute_autocovariance(obj.x_samples_, obj.Sigma);
Gu = obj.compute_autocovariance(obj.u_samples_, obj.Sigma);

G = Gx.*Gu;

W = G + obj.lambda_*M*eye(M);

% Compute value functions.
Vk = zeros(N, M);

in_safe_set = problem.in_safe_set(obj.y_samples_);

cxy = obj.compute_cross_covariance(obj.x_samples_, obj.y_samples_, obj.Sigma);
cuv = obj.compute_cross_covariance(obj.u_samples_, z, obj.Sigma);
beta = cxy.*cuv;
beta = W\beta;
beta = KReach.normalize_beta(beta);

Vk(N, :) = problem.in_target_set(obj.y_samples_);

for k = N-1:-1:2
  Vk(k, :) = in_safe_set.*(Vk(k+1, :)*beta);
end

% Compute probabilities for point.
Pr = zeros(N, mt);

Pr(N, :) = problem.in_target_set(x0);

in_safe_set = problem.in_safe_set(x0);

cxt = obj.compute_cross_covariance(obj.x_samples_, x0, obj.Sigma);
cut = obj.compute_cross_covariance(obj.u_samples_, Ue, obj.Sigma);
beta = cxt.*cut;
beta = W\beta;
beta = obj.normalize_beta(beta);

for k = N-1:-1:1
  Pr(k, :) = in_safe_set.*(Vk(k+1, :)*beta);
end

result = struct;
result.Pr = Pr;

end
