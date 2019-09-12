classdef Tube
% TUBE Tube class.

    properties (Access = private)
        % Sets.
        tube_(:, 1) {mustBeValidSet} = Polyhedron.empty;
    end

    properties (Dependent, Hidden)
        % TUBE The tube.
        tube
    end

    methods
        function obj = Tube(varargin)
            % TUBE Construct an instance of the tube.

            valset = @(arg) mustBeValidSet(arg);

            if nargin
                if nargin == 1
                    % Tube(S)
                    p = inputParser;
                    addOptional(p, 'S', Polyhedron.empty, valset);
                    parse(p, varargin{:});

                    obj.tube_ = p.Results.S;

                elseif nargin == 2
                    % Tube(N, S)
                    p = inputParser;
                    addOptional(p, 'N', 1);
                    addOptional(p, 'S', Polyhedron.empty, valset);
                    parse(p, varargin{:});

                    N = p.Results.N;
                    S = p.Results.S;

                    assert(length(S) == 1);

                    obj.tube_ = repmat(S, N, 1);

                else
                    % Tube(S...)
                    S = [varargin{:}];
                    valset(S);

                    obj.tube_ = S;

                end
            else
                error('Tube is improperly declared. For help, type: help Tube');
            end

        end
    end

    methods
        function tube = get.tube(obj)
            tube = obj.tube_;
        end
        function obj = set.tube(obj, S)
            obj.tube_ = S;
        end
    end

    methods
        function varargout = subsref(obj, s)
            if strcmp(s(1).type, '()')
                [varargout{1:nargout}] = subsref(obj.tube, s);
            else
                [varargout{1:nargout}] = builtin('subsref', obj, s);
            end
        end

        function len = length(obj)
            % LENGTH Returns the length of the tube.
            len = length(obj.tube_);
        end

        function tf = contains(obj, k, varargin)
            if nargin == length(obj.tube_)
                tf = obj.tube_(end).contains(varargin{:});
            else
                tf = obj.tube_(k).contains(varargin{:});
            end
        end

        function v = end(obj, ~, ~)
            v = length(obj);
        end
    end

end

function mustBeValidSet(sets)
% MUSTBEVALIDTUBE Verify tube is valid.

% Sets must be valid MPT classes.
validateattributes(sets, {'Polyhedron', 'Function'}, {'vector'});

% Cannot be empty set.
if ~isempty(sets)
    % for s = sets.'
    %     assert(~s.isEmptySet());
    % end
    sets.forEach(@isEmptySet);
end

end
