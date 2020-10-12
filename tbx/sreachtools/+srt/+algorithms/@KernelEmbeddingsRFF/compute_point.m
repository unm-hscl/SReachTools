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
wx = (1/obj.Sigma)*randn(obj.D, size(sys.X, 1));
wu = (1/obj.Sigma)*randn(obj.D, size(sys.U, 1));
wy = (1/obj.Sigma)*randn(obj.D, size(sys.X, 1));

bx = (2*pi)*rand(obj.D, 1);
bu = (2*pi)*rand(obj.D, 1);
by = (2*pi)*rand(obj.D, 1);

Zx = (sqrt(2)/sqrt(obj.D))*cos(wx*sys.X + bx);
Zu = (sqrt(2)/sqrt(obj.D))*cos(wu*sys.U + bu);
Zy = (sqrt(2)/sqrt(obj.D))*cos(wx*sys.Y + bx);

Z = Zx;

H = Z*Z.';
W = H + obj.Lambda*M*eye(obj.D); %#ok<*MINV>

beta = W\Z;
gamma = beta.'*Zy;

gamma = gamma./sum(gamma, 1);

% Compute value functions.
Vk = zeros(N, M);

switch class(prb)
    case 'srt.problems.FirstHitting'

        Vk(N, :) = prb.TargetTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.TargetTube.contains(k, sys.Y) + ...
                       (prb.ConstraintTube.contains(k, sys.Y) & ...
                        ~prb.TargetTube.contains(k, sys.Y)).* ...
                            (Vk(k+1, :)*gamma);
        end

    case 'srt.problems.TerminalHitting'

        Vk(N, :) = prb.TargetTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.ConstraintTube.contains(k, sys.Y).* ...
                (Vk(k+1, :)*gamma);
        end

    case 'srt.problems.Viability'

        Vk(N, :) = prb.ConstraintTube.contains(N, sys.Y);

        for k = N-1:-1:2
            Vk(k, :) = prb.ConstraintTube.contains(k, sys.Y).* ...
                (Vk(k+1, :)*gamma);
        end

end

% Compute probabilities for point.
Pr = zeros(N, mt);

Zx0 = (sqrt(2)/sqrt(obj.D))*cos(wx*x0 + bx);
Zu0 = (sqrt(2)/sqrt(obj.D))*cos(wu*u0 + bu);

gamma = beta.'*Zx0;

gamma = gamma./sum(gamma, 1);

switch class(prb)
    case 'srt.problems.FirstHitting'

        Pr(N, :) = prb.TargetTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.TargetTube.contains(k, x0) + ...
                       (prb.ConstraintTube.contains(k, x0) & ...
                        ~prb.TargetTube.contains(k, x0)).* ...
                            (Vk(k+1, :)*gamma);
        end

    case 'srt.problems.TerminalHitting'

        Pr(N, :) = prb.TargetTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.ConstraintTube.contains(k, x0).* ...
                (Vk(k+1, :)*gamma);
        end

    case 'srt.problems.Viability'

        Pr(N, :) = prb.ConstraintTube.contains(N, x0);

        for k = N-1:-1:1
            Pr(k, :) = prb.ConstraintTube.contains(k, x0).* ...
                (Vk(k+1, :)*gamma);
        end

end

t_elapsed = toc(t_start);

results = struct;
results.Pr = real(Pr);
results.time = t_elapsed;

end
