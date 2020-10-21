%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Profile SReachOpen algorithms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef ProfileSReachOpen < matlab.perftest.TestCase
% SReachOpen profiles

properties

    Problem

    System

    InitialCondition

end

methods (TestMethodSetup)

    function defineProblem(testCase)

        N = 5;
        K = srt.Tube(N, Polyhedron('lb', [-1; -1], 'ub', [1; 1]));
        T = srt.Tube(N, Polyhedron('lb', [-1; -1], 'ub', [1; 1]));

        testCase.Problem = srt.problems.TerminalHitting( ...
            'ConstraintTube', K, ...
            'TargetTube', T);

    end

    function defineSystem(testCase)

        InputSpace = Polyhedron('lb', 0, 'ub', 0);
        F = eye(2);

        mu = [0; 0];
        sigma = [0.01 0; 0 0.01];

        w = srt.disturbances.Gaussian(mu, sigma);

        testCase.System = srtNDIntegratorModel(2, 0.25, ...
            'F', F, ...
            'w', w, ...
            'InputSpace', InputSpace);

    end

    function defineInitialCondition(testCase)

        testCase.InitialCondition = [0; 0];

    end

end

methods (Test, TestTags = {'NoGurobi'})

    function profileChanceAffine(testCase)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>

        system      = testCase.System;
        x0          = testCase.InitialCondition;

        algorithm   = srt.algorithms.ChanceAffine();

        testCase.startMeasuring();

        SReachPoint(problem, algorithm, system, x0);

        testCase.stopMeasuring();

    end

    function profileChanceAffineUniform(testCase)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>

        system      = testCase.System;
        x0          = testCase.InitialCondition;

        algorithm   = srt.algorithms.ChanceAffineUniform();

        testCase.startMeasuring();

        SReachPoint(problem, algorithm, system, x0);

        testCase.stopMeasuring();

    end

    function profileChanceOpen(testCase)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>

        system      = testCase.System;
        x0          = testCase.InitialCondition;

        algorithm   = srt.algorithms.ChanceOpen();

        testCase.startMeasuring();

        SReachPoint(problem, algorithm, system, x0);

        testCase.stopMeasuring();

    end

    function profileGenzpsOpen(testCase)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>

        system      = testCase.System;
        x0          = testCase.InitialCondition;

        algorithm   = srt.algorithms.GenzpsOpen();

        testCase.startMeasuring();

        SReachPoint(problem, algorithm, system, x0);

        testCase.stopMeasuring();

    end

    function profileKernelEmbeddings(testCase)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>

        A = testCase.System.A;
        B = testCase.System.B;

        s = linspace(-1.1, 1.1, 50);
        X = sampleunif(s, s);

        U = zeros(1, size(X, 2));

        W = 0.01.*randn(size(X));

        Y = A*X + B*U + W;

        args = {'X', X, 'U', U, 'Y', Y};
        samples = srt.systems.SampledSystem(args{:});

        args = {'sigma', 0.1, 'lambda', 1};
        algorithm = srt.algorithms.KernelEmbeddings(args{:});

        x0 = testCase.InitialCondition;
        u0 = 0;

        testCase.startMeasuring();

        SReachPoint(problem, algorithm, samples, x0, u0);

        testCase.stopMeasuring();

    end

    function profileKernelEmbeddingsRFF(testCase)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>

        A = testCase.System.A;
        B = testCase.System.B;

        s = linspace(-1.1, 1.1, 50);
        X = sampleunif(s, s);

        U = zeros(1, size(X, 2));

        W = 0.01.*randn(size(X));

        Y = A*X + B*U + W;

        args = {'X', X, 'U', U, 'Y', Y};
        samples = srt.systems.SampledSystem(args{:});

        args = {'sigma', 0.1, 'lambda', 1, 'D', 5000};
        algorithm = srt.algorithms.KernelEmbeddingsRFF(args{:});

        x0 = testCase.InitialCondition;
        u0 = 0;

        testCase.startMeasuring();

        SReachPoint(problem, algorithm, samples, x0, u0);

        testCase.stopMeasuring();

    end

end

methods (Test)

    function profileParticleOpen(testCase)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>

        system      = testCase.System;
        x0          = testCase.InitialCondition;

        algorithm   = srt.algorithms.ParticleOpen();

        testCase.startMeasuring();

        SReachPoint(problem, algorithm, system, x0);

        testCase.stopMeasuring();

    end

    function profileVoronoiOpen(testCase)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>

        system      = testCase.System;
        x0          = testCase.InitialCondition;

        algorithm   = srt.algorithms.VoronoiOpen();

        testCase.startMeasuring();

        SReachPoint(problem, algorithm, system, x0);

        testCase.stopMeasuring();

    end

end

end
