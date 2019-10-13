classdef Ellipsoid

    properties (Dependent)
        Center

        Sigma
    end

    properties (Dependent, Access = private)
        n_
    end


    properties (Access = private)
        center_

        sigma_
    end

    methods
        function obj = Ellipsoid(center, sigma)
        
            % Input parsing
            p = inputParser();
            p.addRequired('center', @(x) validateattributes(x,...
                {'numeric'}, {'vector','nonempty'}));
            p.addRequired('sigma', @(x) validateattributes(x,...
                {'numeric'}, {'square','nonempty'}));

            p.parse(center, sigma);
            
            center = reshape(center, [], 1);

            obj.center_ = center;
            obj.sigma_ = sigma;

            assert(size(center, 1) == size(sigma, 1));
            assert(issymmetric(sigma));

        end

        function val = get.n_(obj)
            val = length(obj.center_);
        end

        function val = get.Center(obj)
            val = obj.center_;
        end

        function val = get.Sigma(obj)
            val = obj.sigma_;
        end
        
        function val = support(obj, l)

            if isvector(l) && obj.n_ > 1
                reshape(l, [], 1);
            end

            if size(l, 1) ~= obj.n_
                error('Input vector is not of appropriate size.');
            end

            val = zeros(size(l, 2), 1);
            c   = obj.center_;
            sig = obj.sigma_;
            for lv = 1:size(l, 2);
                val(lv) = l(:, lv)' * c + sqrt(l(:, lv)' * sig * l(:, lv));
            end
        end
        
        function newobj = mtimes(obj, F)
        % Override of MATLAB multiplication command
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj - SReachEllipsoid object
        %   F   - Linear transformation matrix for multiplication
        %
        % Outputs:
        % --------
        %   newobj - SReachEllipsoid object (F*obj)
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %
            
            switch [class(obj), class(F)]
                case ['SReachEllipsoid','double']
                    % All ok
                case ['double', 'SReachEllipsoid']
                    % Need to switch the arguments
                    Ftemp = obj;
                    obj = F;
                    F = Ftemp;
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(['Operation *',...
                       ' not defined between *%s, %s'], class(obj), class(F))));
            end
            newobj = srt.Ellipsoid(F * obj.center_, F * obj.sigma_ * F');            
        end
        
        function newobj = plus(obj, v)
        % Override of MATLAB plus command
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj - Ellipsoid object
        %   v   - Deterministic vector to be added to the random vector OR
        %         a Polytope object
        %
        % Outputs:
        % --------
        %   newobj - Ellipsoid obj (obj + v) for deterministic vector/scalar v
        %            Polyhedron obj (obj \oplus v) for polytopic v (overapprox)
        %
        % Notes:
        % ------
        % * For a polytopic v, newobj is an (Polyhedron overapproximation of the
        %   minkowski sum, computed via sampling the support function.
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %
            
            summands_type = [];
            switch [class(obj), class(v)]
                case ['SReachEllipsoid','double']
                    % Check dimensions
                    if ~isequal(size(v), [obj.dim 1]) && ~isequal(size(v),[1 1])
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the SReachEllipsoid and v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'determ_vec_plus_ell';
                case ['double', 'SReachEllipsoid']
                    % Need to switch the arguments
                    vtemp = obj;
                    obj = v;
                    v = vtemp;
                    % Check dimensions
                    if ~isequal(size(v), [obj.dim 1]) && ~isequal(size(v),[1 1]) 
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the random vector and v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'determ_vec_plus_ell';
                case ['SReachEllipsoid','Polyhedron']                
                    % Check dimensions
                    if v.Dim ~= obj.dim
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the SReachEllipsoid and ',...
                            'Polyhedron v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'ell_plus_poly';
                case ['Polyhedron','SReachEllipsoid']                
                    % Need to switch the arguments
                    vtemp = obj;
                    obj = v;
                    v = vtemp;
                    % Check dimensions
                    if v.Dim ~= obj.dim
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the SReachEllipsoid and ',...
                            'Polyhedron v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'ell_plus_poly';
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(['Operation +',...
                       ' not defined between %s and %s'], class(obj),...
                       class(v))));
            end
            switch summands_type
                case 'determ_vec_plus_ell'
                    if isequal(size(v),[1 1])
                        % F is a scalar
                        v = v * ones(obj.dim,1);
                    end
                    newobj = SReachEllipsoid(obj.center_ + v, obj.sigma_);
                case 'ell_plus_poly'
                    new_b = v.support(v.A') + obj.support(v.A');
                    newobj = Polyhedron('H',[v.A new_b]);
                otherwise
                    % Will never come here
            end
        end
        
        function bl = contains(obj, x)
        % Checks if a point (column vector) or a collection of points (matrix of
        % column vectors) is within the ellipsoid
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj         - Ellipsoid object
        %   test_points - Point (column vector) or a N collection of points
        %                 (matrix of column vectors) is within the ellipsoid
        %
        % Outputs:
        % --------
        %   newobj      - Boolean vector Nx1 that describe the containment
        %
        % Notes:
        % ------
        % * Requires CVX for vectorized norm.
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %

            if isvector(x) && obj.n_ > 1
                reshape(x, [], 1);
            end

            if size(x, 1) ~= obj.n_
                error('Input vector is not of appropriate size.');
            end

            % bl = ((x - obj.center_)' / obj.sigma_) * (x - obj.center_) <= 1;
            bl  = zeros(1, size(x, 2));
            c   = obj.center_;
            sig = obj.sigma_;
            for lv = 1:size(x, 2)
                y = x(:, lv);
                bl(lv) = ((y - c)' / sig) * (y - c) <= 1;
            end
        end
    end
end
