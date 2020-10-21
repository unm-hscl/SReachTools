%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Profile KernelEmbeddingsRFF algorithms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef ProfileKernelEmbeddingsRFF < matlab.perftest.TestCase
% KernelEmbeddingsRFF profiles

properties

    Problem

    System
    Samples

    InitialCondition

end

properties (MethodSetupParameter)

    % Number of samples used to build the approximation.
    M = num2cell(100:100:2500);

end

properties (TestParameter)

    % Number of frequency samples used to build the RFF approximation.
    D = num2cell(100:100:1000);

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

    function defineSamples(testCase, M)

        % Generate samples.
        A = testCase.System.A;
        B = testCase.System.B;

        X = -1.1 + 2.2*rand(2, M);

        U = zeros(1, size(X, 2));

        W = 0.01.*randn(size(X));

        Y = A*X + B*U + W;

        args = {'X', X, 'U', U, 'Y', Y};
        testCase.Samples = srt.systems.SampledSystem(args{:});

    end

    function defineInitialCondition(testCase)

        testCase.InitialCondition = [0; 0];

    end

end

methods (Test)

    function profileKernelEmbeddingsRFF(testCase, D)
        % Compute the safety probabilities.
        problem     = testCase.Problem; %#ok<*PROP>
        samples     = testCase.Samples;

        % Define the algorithm.
        args = {'sigma', 0.1, 'lambda', 1, 'D', D};
        algorithm = srt.algorithms.KernelEmbeddingsRFF(args{:});

        % Set the initial condition.
        x0 = testCase.InitialCondition;
        u0 = 0;

        testCase.startMeasuring();

        % Run the test.
        SReachPoint(problem, algorithm, samples, x0, u0);

        testCase.stopMeasuring();

    end

end

end
