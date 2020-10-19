function sys = srtNDIntegratorModel(n, Ts, varargin)

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'n');
addRequired(p, 'Ts');

parse(p, n, Ts, varargin{:});

% anonymous function for getting the necessary matrix internals
facT = @(t, n) t^n / factorial(n);

% initialization
A = eye(n);
B = zeros(n, 1);

% Populate the upper triangle of A and the entries of B
for p = 1:n

    B(p) = facT(Ts, n - p + 1);

    for q = p+1:n
        A(p, q) = facT(Ts, q - p);
    end

end

sys = srt.systems.LtiSystem( ...
    'A', A, ...
    'B', B, ...
    varargin{:} ...
    );

end
