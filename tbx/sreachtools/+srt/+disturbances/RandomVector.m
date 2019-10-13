classdef RandomVector
    properties (Access = protected)
        sample_fun_
    end

    methods
        function obj = RandomVector(varargin)
            p = inputParser();
            addOptional(p, 'sample_fun', @(n) [], ...
                @(x) isa(x, {'function_handle'}));
            parse(p, varargin{:});

            obj.sample_fun_ = p.Results.sample_fun;
        end

        function s = sample(obj, varargin)
            p = inputParser();
            addOptional(p, 'n', 1, @(x) validateattributes(x, {'numeric'}, ...
                {'scalar'}));
            parse(p, varargin{:});

            s = obj.sample_fun_(p.Results.n);
        end

        function rv = concat(obj)
            rv = srt.disturbances.RandomVector();
        end

        function rv = plus(obj, b)
            rv = srt.disturbances.RandomVector(@(n) obj.sample_fun_(n) + b);
        end

        function rv = minus(obj, b)
            rv = srt.disturbances.RandomVector(@(n) obj.sample_fun_(n) - b);
        end

        function rv = times(obj, b)
            rv = srt.disturbances.RandomVector(@(n) b * obj.sample_fun_(n));
        end

        function rv = mtimes(obj, B)
            s = obj.sample();
            if size(B, 2) ~= length(s);
                error(['Incorrect dimensions for matrix multiplication. ', ...
                    'Check that the number of columns in the matrix ', ...
                    'matches the lenght of a sample. To perform ', ...
                    'elementwise multiplication, use ''.*''.']);
            end
            rv = srt.disturbances.RandomVector(@(n) B * obj.sample_fun_(n));
        end
    end

    % methods (Access = protected)
    %     function validatesample(obj)
    %         n = randi([1, 10]);
    %         s = obj.sample(n);
    %         if size(s, 2) ~= n
    %             error(['Sample function is not valid. Sampling functions ', ...
    %                 'must return vectors (nx1).']);
    %         end
            
    %         obj.n_ = size(s, 1);
    %     end
    % end
end