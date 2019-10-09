classdef LagrangianUnder < srt.algorithms.Algorithm
% LAGRANGIANUNDER Lagrangian underapproximation.
    properties (Dependent)
        InputSpace
        
        ProbabilityThreshold

        Method

        VertexFacetAlgorithm
    end

    properties
        input_space_

        prob_threshold_

        method_

        vfalg_
    end

    methods
        function obj = LagrangianUnder(varargin)
            % LAGRANGIANUNDER Construct an instance of the algorithm.

            p = inputParser();
            addParameter(p, 'ProbabilityThreshold', [], ...
                @(x) validateattributes(x, {'numeric'}, {'positive', '<=', 1}));
            addParameter(p, 'InputSpace', [], @(x) isa(x, 'Polyhedron'));
            addParameter(p, 'verbose', 0);
            addParameter(p, 'Method', 'FullEnumeration', ...
                @(x) validatestring(x, {'FullEnumeration', ...
                    'VertexUndersampling'});
            addParameter(p, 'VertexFacetAlgorithm', 'cdd', ...
                @(x) validatestring(x, {'cdd', 'lrs'}));
            parse(p, varargin{:});
            
            if isempty(p.Results.ProbabilityThreshold) || ...
                isempty(p.Results.InputSpace)
                error(['InputSpace and ProbabilityThreshold parameters ', ...
                    'must be specified.']);
            end
            % Call the parent constructor.
            obj = obj@srt.algorithms.Algorithm(varargin{:});

            obj.input_space_ = p.Results.InputSpace;
            obj.prob_threshold_ = p.Results.ProbabilityThreshold;
            obj.vfalg_ = p.Results.VertexFacetAlgorithm;
            obj.method_ = p.Results.Method;

            obj.validateproblem();

        end

        function [R0, varargout] = compute_set(obj, prob, sys, varargin)

            % General process for determining solution
            % 1) Determine the bounded set necessary for the problem
            % 2) Implement recursion
            %   2a) If FullEnumeration, implement recursion exactly as specified
            %   2b) If VertexUndersampling, perform the undersampling
            % 

            time_horizon = prob.TimeHorizon;
            res = struct();

            if strcmp(obj.Method, 'FullEnumeration') && ...
                sys.StateDimension > 3
                
                warning(['You have specified full vertex facet ', ...
                    'enumeration for a system with state dimension > 3. ', ...
                    'Solutions may require exceeding high computing ', ...
                    'resources including time, memory, and processing ', ...
                    'power. Results may also be unstable/irreproducible.']);
                pause(1);
            end

            if sys.StateDimension > 5 && ...
                strcmp(obj.Method, 'VertexUndersampling')

                warning(['You are solving for a higher-dimensional ', ...
                    'problem, with state dimension > 5. Solutions for ', ...
                    'large systems require significant time and ', ...
                    'computational resources.']);
                pause(1);
            end

            % Determine bounded set
            target_set = prob.TargetTube(end);

            % Recursion
            R = prob.TargetTube(end);
            res.ReachTube(time_horizon) = R;
            res.RecursionStepTimer = zeros(1, time_horizon-1);
            full_recursion_timer = tic;
            obj.print_verbose(1, 'Computing recursion.\n')
            for k = time_horizon-1:-1:1
                step_timer = tic;
                if strcmp(obj.Method, 'FullEnumeration')
                    R = full_recurse(k, R, prob.ConstraintTube(k), E, sys);
                else
                    R = vertex_undersample_recurse(k, R, ...
                        prob.ConstraintTube(k), E, sys)
                end
                res.RecursionStepTimer(k) = toc(step_timer);
                res.ReachTube(k) = R;
            end
            res.FullRecursionComputationTime = toc(full_recursion_timer);
            obj.print_verbose(1, ['    Computation time for complete ', ...
                'recursion: %.6f'], res.FullRecursionComputationTime);

            if nargout < 1
                error('Too few output arguments');
            if nargout > 2
                error('Too many output arguments');
            else
                varargout{1} = res;
            end
        end
    end

    methods (Access = private)
        function Rm = full_recurse(k, R, K, E, sys)
            Rm = R - sys.F(k) * E;
            Rm = sys.A(k) \ (Rm + (- sys.B(k) * obj.InputSpace));
            Rm = instersect(Rm, K);
        end

        function Rm = vertex_undersample_recurse(k, R, K, E, sys)
            
        end
    end
    
    methods (Hidden, Access = private)
        function validatedependencies(obj)
            % Needs MPT
            pth = which('Polyhedron')
    
            if isempty(pth)
                error(sprintf(['The Model Parametric toolbox is required ', ...
                    'to use Lagrangian algorithms. You can download it ', ...
                    'at:\n\n', ...
                    '    <a href="https://www.mpt3.org/Main/Installation">https://www.mpt3.org/Main/Installation</a>\n\n']));
            end
        end

        function validateproblem(obj)
            if isempty(obj.InputSpace) || isempty(obj.ProbabilityThreshold)
                error('InputSpace and ProbabilityThreshold must be provided.');
            end
        end
    end
end
