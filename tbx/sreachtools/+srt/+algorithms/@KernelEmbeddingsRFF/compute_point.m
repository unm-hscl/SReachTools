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

% Compute random fourier features.
wx = randn(obj.D, size(sys.X, 1));
% wx = normrnd(0, obj.Sigma^2, obj.D, size(sys.X, 1));
wu = randn(obj.D, size(sys.U, 1));
% wu = normrnd(0, obj.Sigma^2, obj.D, size(sys.U, 1));

b = (2*pi).*rand(obj.D, 1);

Zx = sqrt(2).*cos(wx*sys.X + b);
Zu = sqrt(2).*cos(wu*sys.U + b);

Gx = Zx*Zx.';
Gu = Zu*Zu.';

G = Gx.*Gu;

W = G + obj.lambda_*M*eye(obj.D);

% Compute value functions.
Vk = zeros(N, M);

Zy = sqrt(2).*cos(wx*sys.Y + b);

cxy = Zy;
cuv = Zu;

beta = cxy.*cuv;
beta = (Zx.*Zu).'/W*beta;
beta = obj.normalize_beta(beta);

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

cxt = sqrt(2).*cos(wx*x0 + b);
cut = sqrt(2).*cos(wu*u0 + b);

beta = cxt.*cut;
beta = (Zx.*Zu).'/W*beta;
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
