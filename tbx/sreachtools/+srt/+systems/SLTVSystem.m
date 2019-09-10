classdef SLTVSystem < handle
    properties (SetAccess = protected)
        % System disturbance
        w

        % State space
        X

        % Input space
        U
    end

    properties (Access = protected)
        % State Matrix
        A_

        % Input Matrix
        B_

        % Disturbance Matrix
        F_

        % State dimension
        n_

        % Input dimension
        p_

        % Disturbance dimension
        q_
    end

    methods
        function obj = SLTVSystem(A, B, F, w, U)

            p = inputParser();

            if ~isa(A, 'function_handle')
                validateattributes(A, {'numeric'}, {'square'});
            end

            if ~isa(B, 'function_handle')
                validateattributes(A, {'numeric'}, {'square'});
            end

            if ~isa(F, 'function_handle')
                validateattributes(A, {'numeric'}, {'square'});
            end

            assert(isa(w, 'srt.disturbances.RandomVector'), ...
                ['Disturbances must be of srt.disturbances.RandomVector ', ...
                 'type.']);
            assert(isa(U, 'srt.spaces.Base') || isa(U, 'Polyhedron'), ...
                ['Input space must be of srt.spaces.Base type or ', ...
                 'Polyhedron objects.']);

            obj.A_ = A;
            obj.B_ = B;
            obj.F_ = F;
            obj.w  = w;
            obj.U  = U;

            obj.n_ = size(obj.A(1), 2);
            obj.p_ = size(obj.B(1), 2);
            obj.q_ = size(obj.F(1), 2);

            obj.X = srt.spaces.Rn(obj.n_);
        end

        function xp = onestep(obj, k, x, varargin)
            p = inputParser()
            addRequired(p, 'k', @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
            addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
            addOptional(p, 'u', [], ...
                @(x) validateattributes(x, {'numeric'}, {'vector'}));
            addParameter(p, 'UseMeanValue', false, @(x) islogical(x));
            parse(p, k, x, varargin{:});

            u = p.Results.u;

            Ax = obj.A(k) * x;
            if p.Results.UseMeanValue
                Fw = obj.F(k) * obj.w.mean();
            else
                Fw = obj.F(k) * obj.w.sample();
            end
            Bu = obj.B(k) * u;

            if isempty(Ax) && isempty(Bu) && isempty(Fw)
                xp = [];
            else
                if isempty(Ax) Ax = 0; end
                if isempty(Bu) Bu = 0; end
                if isempty(Fw) Fw = 0; end
                xp = Ax + Bu + Fw;
            end
        end

        function xm = meanstep(obj, k, x, varargin)
            p = inputParser()
            addRequired(p, 'k', @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
            addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
            addOptional(p, 'u', [], ...
                @(x) validateattributes(x, {'numeric'}, {'vector'}));
            parse(p, k, x, varargin{:});

            xm = obj.onestep(k, x, p.Results.u, 'UseMeanValue', true)
        end

        function val = A(obj, k)
            if isa(obj.A_, 'function_handle')
                val = obj.A_(k);
            else
                val = obj.A_;
            end
        end

        function val = B(obj, k)
            if isa(obj.B_, 'function_handle')
                val = obj.B_(k);
            else
                val = obj.B_;
            end
        end

        function val = F(obj, k)
            if isa(obj.F_, 'function_handle')
                val = obj.F_(k);
            else
                val = obj.F_;
            end
        end

        function disp(obj)
            fprintf('  <a href="matlab:helpPopup srt.systems.SLTISystem">SLTVSystem</a> with:\n\n');
            
            if isa(obj.A_, 'function_handle')
                fprintf('    A: [%dx%d time-varying matrix]\n', size(obj.A(1)));
            else
                fprintf('    A: [%dx%d time-invariant matrix]\n', size(obj.A(1)));
            end

            if isa(obj.B_, 'function_handle')
                fprintf('    B: [%dx%d time-varying matrix]\n', size(obj.B(1)));
            else
                fprintf('    B: [%dx%d time-invariant matrix]\n', size(obj.B(1)));
            end

            if isa(obj.F_, 'function_handle')
                fprintf('    F: [%dx%d time-varying matrix]\n', size(obj.F(1)));
            else
                fprintf('    F: [%dx%d time-invariant matrix]\n', size(obj.F(1)));
            end

            fprintf('    w: %s\n', class(obj.w));
            fprintf('    X: %s\n', class(obj.X));
            fprintf('    U: %s\n', class(obj.U));

            fprintf('\n');
        end

        function sys = concat(obj, time_horizon)
        % Create concatenated system

            sys = srt.systems.SLTVSystem( ...
                @(k) obj.getConcatenatedStateMatrix((k-1)*time_horizon + 1: k*time_horizon), ...
                @(k) obj.getConcatenatedInputMatrix((k-1)*time_horizon + 1: k*time_horizon), ...
                @(k) obj.getConcatenatedDisturbanceMatrix((k-1)*time_horizon + 1: k*time_horizon), ...
                obj.w.concat(time_horizon), ...
                obj.concatInputSpace(time_horizon));
        end

        function [Z, H, G] = getConcatenatedMatrices(obj, time_horizon)
        % Get concatenated matrices
        % ============================================================================
        % 
        % Computes the matrices corresponding to the concatentated state vector X.
        %
        % Consider a LtvSystem object with n as the dimension of the state 
        % vector, p as the input dimesnion, and q as the disturbance dimension. 
        % Given a time of interest N, we define a concatenated state vector,
        % an (n*N)-dimensional vector, as
        %           __       __
        %           |   x_1   |
        %           |   x_2   |
        %       X = |   ...   |
        %           | x_{N-1} |          
        %           |  x_{N}  |          
        %           ---     ---
        % where x_t is the state of the system with 1 <= t <= N.  Similarly, one can
        % define concated input and noise vectors U and W (mN-dimensional and
        % pN-dimensional vectors),
        %           __       __         __       __
        %           |   u_0   |         |   w_0   |
        %           |   u_1   |         |   w_1   |
        %       U = |   ...   |,   W  = |   ...   |
        %           | u_{N-2} |         | w_{N-2} |      
        %           | u_{N-1} |         | w_{N-1} |      
        %           ---     ---         ---     ---
        %
        % Given the initial state x_0, we have
        %
        %       X = Z * x_0 + H * U + G * W
        %
        % where Z ((n*N) x n matrix), H ((n*N) x (m*N) matrix), and G  (nN x
        % pN matrix) are appropriate matrices. These matrices (with minor
        % modifications noted below) are given in (3) in 
        %    J. Skaf and S. Boyd, "Design of Affine Controllers via Convex
        %    Optimization", in IEEE Trans. Automatic Control, 2010. 
        %
        % This function computes Z, H, and G.
        %
        % Usage:
        % ------
        %
        % % Compute the concatenated matrices for a double integrator with a time of
        % % interest, 10
        %
        % % Problem parameters
        % time_horizon = 10;
        % T = 0.25;
        % umax = 0.75;
        % dmax = 0.1;
        % % Double integrator system
        % sys = LtvSystem(...
        %     'StateMatrix', [1, T; 0, 1], ...
        %     'InputMatrix', [T^2; T], ...
        %     'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
        %     'DisturbanceMatrix', eye(2), ...
        %     'Disturbance', Polyhedron('lb', -dmax *ones(2,1), 'ub', dmax *ones(2,1)));
        % 
        % % Get the concatenated matrices
        % [Z,H,G] = getConcatMats(sys, time_horizon);
        %
        % =============================================================================
        %
        % [Z,H,G] = getConcatMats(sys, time_horizon)
        %
        % Inputs:
        % -------
        %   sys          - An object of LtvSystem class 
        %   time_horizon - Time horizon (N) with the control provided from 0 to N-1
        %
        % Outputs:
        % --------
        %   Z - Concatenated state matrix
        %   H - Concatenated input matrix
        %   G - Concatenated disturbance matrix
        %
        % Notes:
        % ------
        % * For control-free and/or disturbance-free LTI/LTV systems, H and G are set to
        %   zeros( obj.n_ * time_horizon, 1) as appropriate.
        % * Deviation from Skaf and Boyd's definition,
        %     * Concatenated state is X=[x_1 x_2 ... x_{N}], with the initial state
        %     x_0 EXCLUDED in contrast to Skaf and Boyd's TAC 2010 definitions.
        % * Computes the extended controllability matrix via for loops. (suboptimal way)
        % * This function also serves as a delegatee for input handling
        % 
        % ============================================================================
        %
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %

            % Ensure that time_horizon is a scalar
            assert( isscalar(time_horizon) && time_horizon > 0, ...
                    'SReachTools:invalidArgs', ...
                    'Expected a scalar positive time_horizon');
        
            % Concatenated sizes
            %   dimension x time_horizon
            cn = obj.n_ * time_horizon;
            cp = obj.p_ * time_horizon;
            cq = obj.q_ * time_horizon;
            
            %% Construct Z matrix --- concatenated state matrix
            % Z = [A;A^2;...A^{N}]
            Z = zeros(cn, obj.n_);
            for ith = 1:time_horizon
                Z((ith-1)*obj.n_+1 : ith*obj.n_,:) = ...
                    obj.stateMatrixProduct(0, ith);
            end
        
            %% Construct H matrix --- concatenated input matrix
            if obj.p_ > 0
                H = zeros(cn, cp);
                for ith = 1:time_horizon
                    row_for_H = zeros(obj.n_, cp);
                    for sub_indx = 1:ith
                        input_matrix_temp = obj.B(sub_indx-1);
                        row_for_H(:, ...
                            (sub_indx-1)*obj.p_ + 1: sub_indx*obj.p_) =...
                                obj.stateMatrixProduct(sub_indx, ith) * ...
                                input_matrix_temp;
                    end            
                    H((ith-1)*obj.n_ +1 : ith*obj.n_,:) = row_for_H;
                end
            else
                H = zeros(cn, 0);
            end
        
            %% Construct G matrix --- concatenated disturbance matrix
            if obj.q_ > 0
                G = zeros(cn, cq);
                for t_indx = 1:time_horizon
                    row = zeros(obj.n_, cq);
                    for sub_indx = 1:t_indx
                        dist_matrix_temp = obj.F(sub_indx-1);
                        row(:, (sub_indx-1)*obj.q_ + 1: sub_indx*obj.q_) = ...
                            obj.stateMatrixProduct(sub_indx, t_indx) * ...
                            dist_matrix_temp;
                    end            
                    G((t_indx-1)*obj.n_ +1 : t_indx*obj.n_,:) = row;
                end
            else
                G = zeros(cn, 0);
            end
        end

        function [cA, cb] =  getConcatenatedInputSpace(obj, time_horizon)
        % Get half space representation of the concatenated (polytopic) input space 
        % for the given time horizon
        % ============================================================================
        % 
        % Computes the input_space^{time_horizon} corresponding to a given set, which
        % is the set of admissible open-loop control polices. This function computes the
        % half-space representation of the cartesian products of polytopic input spaces.
        %
        % Usage:
        % ------
        %
        % % Compute the (matrix form) set of admissible open-loop control policies given
        % % a LtvSystem and a time horizon
        %
        % sys = LtvSystem(...
        %     'StateMatrix', eye(2), ...
        %     'InputMatrix', ones(2,1), ...
        %     'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
        % time_horizon = 10;
        % [cA, cb] = ...
        %                                             getConcatenatedInputSpace(sys, ...
        %                                                                 time_horizon);
        % 
        % ============================================================================
        %
        % [cA, cb] =...
        %                                              getConcatenatedInputSpace(sys, ...
        %                                                                  time_horizon)
        % 
        % Inputs:
        % -------
        %   sys                  - An object of LtvSystem class 
        %   time_horizon         - Time horizon
        %
        % Outputs:
        % --------
        %   cA, cb 
        %                        - Concatenated input space (Halfspace representation)
        %
        % =============================================================================
        %
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %

            %% Input handling
            % Ensure that the system has a non-empty input space
            if obj.U.isEmptySet
                throwAsCaller(SrtInvalidArgsError(['Expected a non-empty ', ...
                'polyhedral input space']));
            end

            % Ensure that time horizon is a scalar and positive
            if ~isscalar(time_horizon) || time_horizon <= 0
                throwAsCaller(SrtInvalidArgsError(['Expected a scalar ', ...
                'positive time_horizon']));
            end

            %% Construction of the concatenated input space (input_space^{time_horizon})
            cA = kron(eye(time_horizon), obj.U.A);
            cb = kron(ones(time_horizon,1), obj.U.b);
        end

        function getConcatenatedStateMatrix(obj, t_start, t_end)
            %% Construct Z matrix --- concatenated state matrix
            % Z = [A;A^2;...A^{N}]
            Z = zeros(cn, obj.n_);
            for ith = t_start:t_end
                Z((ith-1)*obj.n_+1 : ith*obj.n_,:) = ...
                    obj.stateMatrixProduct(0, ith);
            end
        end

        function getConcatenatedInputMatrix(obj, t_start, t_end)
            %% Construct H matrix --- concatenated input matrix
            if obj.p_ > 0
                H = zeros(cn, cp);
                for ith = t_start:t_end
                    row_for_H = zeros(obj.n_, cp);
                    for sub_indx = 1:ith
                        input_matrix_temp = obj.B(sub_indx-1);
                        row_for_H(:, ...
                            (sub_indx-1)*obj.p_ + 1: sub_indx*obj.p_) =...
                                obj.stateMatrixProduct(sub_indx, ith) * ...
                                input_matrix_temp;
                    end            
                    H((ith-1)*obj.n_ +1 : ith*obj.n_,:) = row_for_H;
                end
            else
                H = zeros(cn, 0);
            end
        end

        function getConcatenatedDisturbanceMatrix(obj, t_start, t_end)
            %% Construct G matrix --- concatenated disturbance matrix
            if obj.q_ > 0
                G = zeros(cn, cq);
                for t_indx = t_start:t_end
                    row = zeros(obj.n_, cq);
                    for sub_indx = 1:t_indx
                        dist_matrix_temp = obj.F(sub_indx-1);
                        row(:, (sub_indx-1)*obj.q_ + 1: sub_indx*obj.q_) = ...
                            obj.stateMatrixProduct(sub_indx, t_indx) * ...
                            dist_matrix_temp;
                    end            
                    G((t_indx-1)*obj.n_ +1 : t_indx*obj.n_,:) = row;
                end
            else
                G = zeros(cn, 0);
            end
        end
    end

    methods (Access = private)
        function P = stateMatrixProduct(obj, t, tau)
        % Compute product of state matrices between times [t, tau]:
        % 
        %   P = A(t) * A(t+1) * A(t+2) * ... * A(tau)
        % 
        % If t >= tau: P is an identity matrix
        % 

            P = eye(obj.n_);

            if t<tau
                for lv = t:tau-1
                    P = P * obj.A(lv);
                end    
            end
        end
    end
end