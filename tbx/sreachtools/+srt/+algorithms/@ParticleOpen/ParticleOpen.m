classdef ParticleOpen < srt.algorithms.Algorithm
% PARTICLEOPEN Particle control.

    properties (Access = private)

        % Number of particles to be used for approximation
        n_particles (1, 1) double {mustBePositive} = 100

        % BigM notation requires a large value
        bigM (1, 1) double {mustBePositive} = 100

        % Maximum number of particles allowed for tractable computation
        max_particles (1, 1) double {mustBePositive} = 200

    end

    methods
        function obj = ParticleOpen(varargin)
            % PARTICLEOPEN Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            p = inputParser;
            p.KeepUnmatched = true;

            addParameter(p, 'n_particles', 100);
            addParameter(p, 'bigM', 100);
            addParameter(p, 'max_particles', 200);

            parse(p, varargin{:});

            obj.n_particles = p.Results.n_particles;
            obj.bigM = p.Results.bigM;
            obj.max_particles = p.Results.max_particles;

        end
    end

end
