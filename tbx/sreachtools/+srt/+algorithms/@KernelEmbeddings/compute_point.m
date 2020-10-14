function results = compute_point(obj, prob, sys, x0, u0, varargin)

p = inputParser;
addRequired(p, 'prob', @obj.validateproblem);
addRequired(p, 'sys', @obj.validatesystem);
addRequired(p, 'x0');
parse(p, prob, sys, x0, varargin{:});

import srt.*

% Constants
N = prob.TimeHorizon;

mt = size(x0, 2);

M = sys.length;

t_start = tic;

% Compute weight matrix.
obj.print_verbose(1, 'Computing weight matrix.\n');

Gx = obj.compute_autocovariance(sys.X, obj.Sigma);
Gu = obj.compute_autocovariance(sys.U, obj.Sigma);

G = Gx.*Gu;

W = inv(G + obj.lambda_*M*eye(M));

% Compute value functions.
obj.print_verbose(1, 'Computing value functions.\n');

Vk = zeros(N, M);

cxy = obj.compute_cross_covariance(sys.X, sys.Y, obj.Sigma);
cuv = obj.compute_cross_covariance(sys.U, sys.U, obj.Sigma);

beta = cxy.*cuv;
beta = W*beta; %#ok<*MINV>
beta = algorithms.KernelEmbeddings.normalize_beta(beta);

switch class(prob)
    case 'srt.problems.FirstHitting'

        Vk(N, :) = prob.TargetTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prob.TargetTube.contains(k, sys.Y) + ...
                       (prob.ConstraintTube.contains(k, sys.Y) & ...
                        ~prob.TargetTube.contains(k, sys.Y)).*(Vk(k+1, :)*beta);
        end

    case 'srt.problems.TerminalHitting'

        Vk(N, :) = prob.TargetTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prob.ConstraintTube.contains(k, sys.Y).*(Vk(k+1, :)*beta);
        end

    case 'srt.problems.Viability'

        Vk(N, :) = prob.ConstraintTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prob.ConstraintTube.contains(k, sys.Y).*(Vk(k+1, :)*beta);
        end

end

% Compute probabilities for point.
obj.print_verbose(1, 'Computing safety probabilities.\n');

Pr = zeros(N, mt);

cxt = obj.compute_cross_covariance(sys.X, x0, obj.Sigma);
cut = obj.compute_cross_covariance(sys.U, u0, obj.Sigma);

beta = cxt.*cut;
beta = W*beta;
beta = obj.normalize_beta(beta);

switch class(prob)
    case 'srt.problems.FirstHitting'

        Pr(N, :) = prob.TargetTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prob.TargetTube.contains(k, x0) + ...
                       (prob.ConstraintTube.contains(k, x0) & ...
                        ~prob.TargetTube.contains(k, x0)).*(Vk(k+1, :)*beta);
        end

    case 'srt.problems.TerminalHitting'

        Pr(N, :) = prob.TargetTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prob.ConstraintTube.contains(k, x0).*(Vk(k+1, :)*beta);
        end

    case 'srt.problems.Viability'

        Pr(N, :) = prob.ConstraintTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prob.ConstraintTube.contains(k, x0).*(Vk(k+1, :)*beta);
        end

end

t_elapsed = toc(t_start);

results = struct;
results.Pr = Pr;
results.time = t_elapsed;

end
