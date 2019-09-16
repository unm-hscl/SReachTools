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

t_start = tic;

% Compute random fourier features.
wx = (1/obj.Sigma^2).*randn(obj.D, size(sys.X, 1));
wu = (1/obj.Sigma^2).*randn(obj.D, size(sys.U, 1));

Zx = exp(1i*wx*sys.X);
Zu = exp(1i*wu*sys.U);
Zy = exp(1i*wx*sys.Y);

Z = Zx.*Zu;

W = Z*Z' + obj.Lambda*eye(obj.D);

beta = W\Z;

% Compute value functions.
Vk = zeros(N, M);

switch class(prb)
    case 'srt.problems.FirstHitting'

        Vk(N, :) = prb.TargetTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.TargetTube.contains(k, sys.Y) + ...
                       (prb.ConstraintTube.contains(k, sys.Y) & ...
                        ~prb.TargetTube.contains(k, sys.Y)).* ...
                            (Vk(k+1, :)*beta'*(Zy.*Zu));
        end

    case 'srt.problems.TerminalHitting'

        Vk(N, :) = prb.TargetTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.ConstraintTube.contains(k, sys.Y).* ...
                (Vk(k+1, :)*beta'*(Zy.*Zu));
        end

    case 'srt.problems.Viability'

        Vk(N, :) = prb.ConstraintTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.ConstraintTube.contains(k, sys.Y).* ...
                (Vk(k+1, :)*beta'*(Zy.*Zu));
        end

end

% Compute probabilities for point.
Pr = zeros(N, mt);

Zx0 = exp(1i*wx*x0);
Zu0 = exp(1i*wu*u0);

switch class(prb)
    case 'srt.problems.FirstHitting'

        Pr(N, :) = prb.TargetTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.TargetTube.contains(k, x0) + ...
                       (prb.ConstraintTube.contains(k, x0) & ...
                        ~prb.TargetTube.contains(k, x0)).* ...
                            (Vk(k+1, :)*beta'*(Zx0.*Zu0));
        end

    case 'srt.problems.TerminalHitting'

        Pr(N, :) = prb.TargetTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.ConstraintTube.contains(k, x0).* ...
                (Vk(k+1, :)*beta'*(Zx0.*Zu0));
        end

    case 'srt.problems.Viability'

        Pr(N, :) = prb.ConstraintTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.ConstraintTube.contains(k, x0).* ...
                (Vk(k+1, :)*beta'*(Zx0.*Zu0));
        end

end

t_elapsed = toc(t_start);

results = struct;
results.Pr = real(Pr);
results.time = t_elapsed;

end
