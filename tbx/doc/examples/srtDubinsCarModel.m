function sys = srtDubinsCarModel(N, Ts, varargin)

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'N');
addOptional(p, 'Ts', 0.1);

addParameter(p, 'InitialHeading', pi/10);

addParameter(p, 'DisturbanceType', 'affine', ...
    @(arg) any(strcmpi(arg, {'affine', 'velocity'})));

% addParameter(p, 'StateSpace', Polyhedron('lb', 0, 'ub', 10));
addParameter(p, 'InputSpace', Polyhedron('lb', 0, 'ub', 10));

addParameter(p, 'Disturbance', ...
    srt.disturbances.Gaussian(zeros(2, 1), 0.001 * eye(2)));

addParameter(p, 'DisturbanceMatrix', eye(2));

parse(p, N, Ts, varargin{:});

InitialHeading = p.Results.InitialHeading;

% StateSpace = p.Results.StateSpace;
InputSpace = p.Results.InputSpace;
Disturbance = p.Results.Disturbance;
DisturbanceMatrix = p.Results.DisturbanceMatrix;

% Known turning rate sequence
omega = pi/N/Ts;
TurningRate = omega*ones(N, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HeadingVector = InitialHeading + Ts * cumsum([0; TurningRate]);

TimeVaryingMatrix = @(t) Ts * [cos(HeadingVector(t + 1)); ...
                               sin(HeadingVector(t + 1))];

if strcmpi(p.Results.DisturbanceType, 'affine')

    sys = srt.systems.LtvSystem( ...
        'A', @(t) eye(2), ...
        'B', TimeVaryingMatrix, ...
        'InputSpace', InputSpace, ...
        'F', DisturbanceMatrix, ...
        'w', Disturbance);

else

    sys = srt.systems.LtvSystem( ...
        'A', @(t) eye(2), ...
        'B', TimeVaryingMatrix, ...
        'InputSpace', InputSpace, ...
        'F', TimeVaryingMatrix, ...
        'w', Disturbance);

end

end
