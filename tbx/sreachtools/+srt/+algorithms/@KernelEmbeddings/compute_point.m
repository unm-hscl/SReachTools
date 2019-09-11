function results = compute_point(obj, prb, sys, x0, u0, varargin)

p = inputParser;
addRequired(p, 'prb', @obj.validateproblem);
addRequired(p, 'sys', @obj.validatesystem);
addRequired(p, 'x0');
parse(p, prb, sys, x0, varargin{:});

% Constants
N = prb.TimeHorizon;

mt = size(x0, 2);

M = sys.length;

% Compute weight matrix.
Gx = obj.compute_autocovariance(sys.X, obj.Sigma);
Gu = obj.compute_autocovariance(sys.U, obj.Sigma);

G = Gx.*Gu;

W = G + obj.lambda_*M*eye(M);

% Compute value functions.
Vk = zeros(N, M);

z = zeros(size(sys.U));

cxy = obj.compute_cross_covariance(sys.X, sys.Y, obj.Sigma);
cuv = obj.compute_cross_covariance(sys.U, z, obj.Sigma);
beta = cxy.*cuv;
beta = W\beta;
beta = srt.algorithms.KernelEmbeddings.normalize_beta(beta);

Vk(N, :) = prb.contains(N, sys.Y);

for k = N-1:-1:2
  Vk(k, :) = prb.contains(k, sys.Y).*(Vk(k+1, :)*beta);
end

% Compute probabilities for point.
Pr = zeros(N, mt);

Pr(N, :) = prb.contains(N, x0);

cxt = obj.compute_cross_covariance(sys.X, x0, obj.Sigma);
cut = obj.compute_cross_covariance(sys.U, u0, obj.Sigma);
beta = cxt.*cut;
beta = W\beta;
beta = obj.normalize_beta(beta);

for k = N-1:-1:1
  Pr(k, :) = prb.contains(k, x0).*(Vk(k+1, :)*beta);
end

results = struct;
results.Pr = Pr;

end
