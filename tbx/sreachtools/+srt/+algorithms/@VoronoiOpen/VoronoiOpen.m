classdef VoronoiOpen < srt.algorithms.Algorithm
% VORONOIOPEN Voronoi partition control.

    properties (Access = private)

        % Risk of the probabilistic overapproximation bound failing
        failure_risk (1, 1) double {mustBePositive} = 1E-4

        % Maximum overapproximation error (probabilistically) tolerable
        max_overapprox_err (1, 1) double {mustBePositive} = 1E-2

        % Number of kmeans cluster points/ Voronoi centers
        n_kmeans (1, 1) double {mustBePositive} = 30

        % Maximum number of particles allowed for tractable computation
        max_particles (1, 1) double {mustBePositive} = 1E5

        % BigM notation requires a large value
        bigM (1, 1) double {mustBePositive} = 100

    end

    methods
        function obj = VoronoiOpen(varargin)
            % VORONOIOPEN Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            p = inputParser;
            p.KeepUnmatched = true;

            addParameter(p, 'failure_risk', 1E-4);
            addParameter(p, 'max_overapprox_err', 1E-2);
            addParameter(p, 'n_kmeans', 30);
            addParameter(p, 'max_particles', 1E5);
            addParameter(p, 'bigM', 100);

            parse(p, varargin{:});

            obj.failure_risk = p.Results.failure_risk;
            obj.max_overapprox_err = p.Results.max_overapprox_err;
            obj.n_kmeans = p.Results.n_kmeans;
            obj.max_particles = p.Results.max_particles;
            obj.bigM = p.Results.bigM;


        end
    end

end
