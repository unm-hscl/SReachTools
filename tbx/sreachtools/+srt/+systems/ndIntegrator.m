function sys = ndIntegrator(n, T, U, varargin)
    p = inputParser();
    addRequired(p, 'n', @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'integer', 'positive'}));
    addRequired(p, 'T', @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'positive'}));
    addRequired(p, 'U', @(x) isa(x, 'srt.spaces.Base'));
    addOptional(p, 'F', [], @(x) validateattributes(x, {'numeric'}, ...
        {'nrows', n}));
    addOptional(p, 'w', srt.disturbances.Empty(), ...
        @(x) isa(x, 'srt.disturbances.RandomVector'));

    parse(p, n, T, U, varargin{:});

    F = p.Results.F;
    w = p.Results.w;

    % anonymous function for getting the necessary matrix internals
    facT = @(t, n) t^n / factorial(n);

    % initialization
    A = eye(n);
    B = zeros(n, 1);
    
    % Populate the upper triangle of A and the entries of B
    for i = 1:n
        B(i) = facT(T, n-i+1);
        for j = i+1:n
            A(i, j) = facT(T, j-i);
        end
    end

    sys = srt.systems.SLTISystem(A, B, F, w, U);
end