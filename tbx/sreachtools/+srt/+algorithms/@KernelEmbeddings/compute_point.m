function results = compute_point(obj, prb, sys, x0, u0, varargin)

p = inputParser;
addRequired(p, 'prb', @obj.validateproblem);
addRequired(p, 'sys', @obj.validatesystem);
addRequired(p, 'x0');
parse(p, prb, sys, x0, varargin{:});

import srt.*

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

cxy = obj.compute_cross_covariance(sys.X, sys.Y, obj.Sigma);
cuv = obj.compute_cross_covariance(sys.U, sys.U, obj.Sigma);
beta = cxy.*cuv;
beta = W\beta;
beta = algorithms.KernelEmbeddings.normalize_beta(beta);

switch class(prb)
    case 'srt.problems.FirstHitting'

        Vk(N, :) = prb.TargetTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.TargetTube.contains(k, sys.Y) + ...
                       (prb.ConstraintTube.contains(k, sys.Y) & ...
                        ~prb.TargetTube.contains(k, sys.Y)).*(Vk(k+1, :)*beta);
        end

    case 'srt.problems.TerminalHitting'

        Vk(N, :) = prb.TargetTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.ConstraintTube.contains(k, sys.Y).*(Vk(k+1, :)*beta);
        end

    case 'srt.problems.Viability'

        Vk(N, :) = prb.ConstraintTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.ConstraintTube.contains(k, sys.Y).*(Vk(k+1, :)*beta);
        end

end

% Compute probabilities for point.
Pr = zeros(N, mt);

cxt = obj.compute_cross_covariance(sys.X, x0, obj.Sigma);
cut = obj.compute_cross_covariance(sys.U, u0, obj.Sigma);
beta = cxt.*cut;
beta = W\beta;
beta = obj.normalize_beta(beta);

switch class(prb)
    case 'srt.problems.FirstHitting'

        Pr(N, :) = prb.TargetTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.TargetTube.contains(k, x0) + ...
                       (prb.ConstraintTube.contains(k, x0) & ...
                        ~prb.TargetTube.contains(k, x0)).*(Vk(k+1, :)*beta);
        end

    case 'srt.problems.TerminalHitting'

        Pr(N, :) = prb.TargetTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.ConstraintTube.contains(k, x0).*(Vk(k+1, :)*beta);
        end

    case 'srt.problems.Viability'

        Pr(N, :) = prb.ConstraintTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.ConstraintTube.contains(k, x0).*(Vk(k+1, :)*beta);
        end

end

results = struct;
results.Pr = Pr;

end
