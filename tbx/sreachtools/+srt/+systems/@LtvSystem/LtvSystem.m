classdef LtvSystem < srt.systems.StochasticSystem
    properties (Dependent)
        % STATESPACE State Space of the system
        StateSpace

        % DISTURBANCE Disturbance object for the system
        Disturbance

        % STATEDIMENSION Dimension of the state
        StateDimension

        % INPUTDIMENSION Dimension of the input
        InputDimension

        % DISTURBANCEDIMENSION Dimension of the disturbance
        DisturbanceDimension
    end

    properties (Access = protected)
        % A_ State Matrix
        A_

        % B_ Input Matrix
        B_

        % F_ Disturbance Matrix
        F_

        % W_ Disturbance
        w_

        % X_ State Space
        X_

        % N_ State dimension
        n_

        % P_ Input dimension
        p_

        % Q_ Disturbance dimension
        q_
    end

    methods
        function obj = LtvSystem(varargin)
            obj = obj@srt.systems.StochasticSystem();

            p = inputParser();
            addOptional(p, 'A', [], @(x) validateattributes(x, ...
                {'function_handle', 'numeric'}, {'square'}));
            addOptional(p, 'B', [], ...
                @(x) isa(x, 'function_handle') || isa(x, 'numeric'));
            addOptional(p, 'F', [], ...
                @(x) obj.validateDisturbanceMatrix(x));
            addOptional(p, 'w', srt.disturbances.RandomVector(), ...
                @(x) isa(x, 'srt.disturbances.RandomVector'));
            parse(p, varargin{:});

            obj.A_ = p.Results.A;
            obj.B_ = p.Results.B;
            if ischar(p.Results.F)
                obj.F_ = p.Results.B;
            else
                obj.F_ = p.Results.F;
            end
            obj.w_ = p.Results.w;

            obj.X_ = srt.spaces.Rn(obj.n_);

            obj.n_ = size(obj.A(1), 2);

            if size(obj.B(1), 1) > 0 && size(obj.B(1), 1) ~= obj.n_
                error(['Input matrix must be empty or have same number ', ...
                    'of rows as the state matrix.']);
            end
            obj.p_ = size(obj.B(1), 2);

            if size(obj.F(1), 1) > 0 && size(obj.F(1), 1) ~= obj.n_
                error(['Disturbance matrix must be empty or have same ', ...
                    'number of rows as the state matrix.']);
            end
            if size(obj.F(1), 2) ~= length(obj.w_.sample())
                error(['Disturbance matrix and disturbance dimensions do ', ...
                    'not match, i.e. size(F, 2) ~= length(w.sample()).']);
            end
            obj.q_ = size(obj.F(1), 2);

        end

        function val = get.StateSpace(obj)
            val = obj.X_;
        end

        function val = get.Disturbance(obj)
            val = obj.w_;
        end

        function val = get.StateDimension(obj)
            val = obj.n_;
        end

        function val = get.InputDimension(obj)
            val = obj.p_;
        end

        function val = get.DisturbanceDimension(obj)
            val = obj.q_;
        end

        function xp = onestep(obj, k, x, varargin)
            p = inputParser();
            addRequired(p, 'k', @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
            addRequired(p, 'x', @(x) validateattributes(x, {'numeric'}, ...
                {'vector', 'numel', obj.n_}));
            addOptional(p, 'u', [], ...
                @(x) validateattributes(x, {'numeric'}, ...
                    {'vector', 'numel', obj.p_}));
            addOptional(p, 'w', [], ...
                @(x) validateattributes(x, {'numeric'}, ...
                    {'vector', 'numel', obj.q_}));
            parse(p, k, x, varargin{:});

            x = reshape(x, [], 1);
            Ax = obj.A(k) * x;

            if obj.p_ > 0 && any(strcmp('u', p.UsingDefaults))
                error(['No input provided to LtvSystem with nonzero ', ...
                    'input dimension']);
            end
            u = reshape(p.Results.u, [], 1);
            Bu = obj.B(k) * u;

            if obj.q_ > 0 && any(strcmp('w', p.UsingDefaults))
                Fw = obj.F(k) * obj.w_;
            else
                Fw = obj.F(k) * p.Results.w;
            end

            if isempty(Ax) && isempty(Bu) && isempty(Fw)
                xp = [];
            else
                if isempty(Ax), Ax = 0; end
                if isempty(Bu), Bu = 0; end
                if isempty(Fw), Fw = 0; end
                xp = Ax + Bu + Fw;
            end
        end

        function val = StateMatrix(obj, k)
            val = obj.A(k);
        end

        function val = A(obj, k)
            if isa(obj.A_, 'function_handle')
                val = obj.A_(k);
            else
                val = obj.A_;
            end
        end

        function val = InputMatrix(obj, k)
            val = obj.B(k);
        end

        function val = B(obj, k)
            if isa(obj.B_, 'function_handle')
                val = obj.B_(k);
            else
                val = obj.B_;
            end
        end

        function val = DisturbanceMatrix(obj, k)
            val = obj.F(k);
        end

        function val = F(obj, k)
            if isa(obj.F_, 'function_handle')
                val = obj.F_(k);
            else
                val = obj.F_;
            end
        end

        function sys = concat(obj, time_horizon)
        % Create concatenated system

            sys = srt.systems.SLTVSystem( ...
                @(k) obj.concatenateStateMatrix((k-1)*time_horizon + 1: k*time_horizon), ...
                @(k) obj.concatenateInputMatrix((k-1)*time_horizon + 1: k*time_horizon), ...
                @(k) obj.concatenateDisturbanceMatrix((k-1)*time_horizon + 1: k*time_horizon), ...
                obj.w.concat(time_horizon));
        end

        function concatenateStateMatrix(obj, t_start, t_end)
            %% Construct Z matrix --- concatenated state matrix
            % Z = [A;A^2;...A^{N}]
            Z = zeros(cn, obj.n_);
            for ith = t_start:t_end
                Z((ith-1)*obj.n_+1 : ith*obj.n_,:) = ...
                    obj.stateMatrixProduct(0, ith);
            end
        end

        function concatenateInputMatrix(obj, t_start, t_end)
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

        function concatenateDisturbanceMatrix(obj, t_start, t_end)
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

    methods (Static)
        function valid = validateDisturbanceMatrix(F)
            if isa(F, 'function_hanlde') || ...
                (isa(F, 'numeric') && length(size(F)) == 2) || ...
                (isa(F, 'char') && any(strcmp(F, {'InputMatrix', 'B'})))

                valid = true;
            else
                error(sprintf(['Expected input to be one of these types:\n\n', ...
                    'function_handle, double, single, uint8, uint16, ', ...
                    'uint32, uint64, int8, int16, int32, int64\n\n', ...
                    'Or should be a character array matching either:\n\n', ...
                    '''InputMatrix'', ''B''']));
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