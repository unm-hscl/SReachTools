classdef LagrangianUnder < srt.algorithms.Algorithm
% LAGRANGIANUNDER Lagrangian underapproximation.
    properties (Dependent)
        Method

        VertexFacetAlgorithm

        TemplatePolytope
    end

    properties (Access = private)
        method_

        vfalg_

        template_poly_
    end

    methods
        function obj = LagrangianUnder(varargin)
            % LAGRANGIANUNDER Construct an instance of the algorithm.

            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            p = inputParser();
            addParameter(p, 'Method', 'FullEnumeration', ...
                @(x) validatestring(x, {'FullEnumeration', ...
                    'VertexUndersampling'}));
            addParameter(p, 'VertexFacetAlgorithm', 'cdd', ...
                @(x) validatestring(x, {'cdd', 'lrs'}));
            addParameter(p, 'TemplatePolytope', Polyhedron(), ...
                @(x) isa(x, 'Polyhedron') || isa(x, 'srt.Ellipsoid'));
            addParameter(p, 'verbose', 0);
            parse(p, varargin{:});


            obj.vfalg_ = p.Results.VertexFacetAlgorithm;
            obj.method_ = p.Results.Method;

            obj.template_poly_ = p.Results.TemplatePolytope;

        end

        function val = get.Method(obj)
            val = obj.method_;
        end

        function val = get.VertexFacetAlgorithm(obj)
            val = obj.vfalg_;
        end

        function val = get.TemplatePolytope(obj)
            val = obj.template_poly_;
        end

        function res = compute_set(obj, prob, sys, lev, U)

            % General process for determining solution
            % 1) Determine the bounded set necessary for the problem
            % 2) Implement recursion
            %   2a) If FullEnumeration, implement recursion exactly as specified
            %   2b) If VertexUndersampling, perform the undersampling
            %

            p = inputParser();
            addRequired(p, 'prob', @(x) obj.validateproblem(x));
            addRequired(p, 'sys', @(x) obj.validatesystem(x));
            addRequired(p, 'lev', @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'nonnegative', '<=', 1}));
            addRequired(p, 'U', @(x) isa(x, 'Polyhedron'));
            parse(p, prob, sys, lev, U);

            time_horizon = prob.TimeHorizon;
            res = struct();

            if strcmp(obj.Method, 'FullEnumeration') && ...
                sys.StateDimension > 3

                warning(['You have specified full vertex facet ', ...
                    'enumeration for a system with state dimension > 3. ', ...
                    'Solutions may require exceedingly high computing ', ...
                    'resources including time, memory, and processing ', ...
                    'power. Results may also be unstable/irreproducible.']);
            end

            if sys.StateDimension > 5 && ...
                strcmp(obj.Method, 'VertexUndersampling')

                warning(['You are solving for a higher-dimensional ', ...
                    'problem, with state dimension > 5. Solutions for ', ...
                    'large systems require significant time and ', ...
                    'computational resources.']);
            end

            % Determine bounded set
            % E = Polyhedron();
            E = obj.compute_bounded_set(time_horizon, lev, sys.Disturbance);

            % Recursion
            R = prob.TargetTube(end);
            res.ReachTube(time_horizon) = R;
            res.RecursionStepTimer = zeros(1, time_horizon-1);
            full_recursion_timer = tic;
            obj.print_verbose(1, 'Computing recursion.\n')
            for k = time_horizon-1:-1:1
                step_timer = tic;
                if strcmp(obj.Method, 'FullEnumeration')
                    R = obj.full_recursion(k, R, prob.ConstraintTube(k), ...
                        E, U, sys);
                else
                    R = obj.vertex_undersample_recurse(k, R, ...
                        prob.ConstraintTube(k), E, U, sys)
                end
                res.RecursionStepTimer(k) = toc(step_timer);
                res.ReachTube(k) = R;
                obj.print_verbose(2, ['    Computataion time for ', ...
                    'k = %d: %.5f\n'], res.RecursionStepTimer);
            end
            res.FullRecursionComputationTime = toc(full_recursion_timer);
            obj.print_verbose(1, ['    Computation time for complete ', ...
                'recursion: %.5f\n'], res.FullRecursionComputationTime);
        end
    end

    methods (Access = private)
        function Rm = full_recursion(obj, k, R, K, E, U, sys)
            if isa(E, 'srt.Ellipsoid')
                Rm = Polyhedron('H', [R.A, R.b - E.support(R.A')]);
            else
                if isEmptySet(E)
                    Rm = R;
                else
                    Rm = R - sys.F(k) * E;
                end
            end

            Rm = inv(sys.A(k)) * (Rm + (- sys.B(k) * U));
            Rm = intersect(Rm, K);
        end

        function Rm = vertex_undersample_recursion(obj, k, R, K, E, U, sys)

        end
    end

    methods (Hidden, Access = private)
        function E = compute_bounded_set(obj, N, lev, dv)
            % Probability level required of the bounded set
            gam = lev^(1/N);

            if isEmptySet(obj.TemplatePolytope) && ...
                isa(dv, 'srt.disturbances.Gaussian')
                % Get ellipsoid through chi-squared inversion
                R2 = chi2inv(gam, length(dv.Mean));

                E = srt.Ellipsoid(dv.Mean, dv.Sigma * R2);
            else
                % Empirical estimate using bisection

                a = 0;
                b = 1;

            end
        end

        function validatedependencies(obj)
            % Needs MPT
            pth = which('Polyhedron');

            if isempty(pth)
                error(sprintf(['The Model Parametric toolbox is required ', ...
                    'to use Lagrangian algorithms. You can download it ', ...
                    'at:\n\n', ...
                    '    <a href="https://www.mpt3.org/Main/Installation">https://www.mpt3.org/Main/Installation</a>\n\n']));
            end
        end

        function validateproblem(obj, prob)
            if ~isa(prob, 'srt.problems.TerminalHitting')
                error('Problem must be TerminalHitting.');
            end

            if ~isPolyhedronTube(prob.TargetTube) || ...
                ~isPolyhedronTube(prob.ConstraintTube)
                error('Target and Constraint tubes must be polyhedral.');
            end
        end

        function validatesystem(obj, sys)
            if ~isa(sys, 'srt.systems.LtvSystem')
                error('Lagrangian methods require LtvSystems or LtiSystems.');
            end
        end
    end
end
