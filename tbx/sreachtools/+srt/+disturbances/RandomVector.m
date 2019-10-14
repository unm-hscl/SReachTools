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

        function rv = plus(A, B)
            if isa(A, 'srt.disturbances.RandomVector')
                rv = srt.disturbances.RandomVector(@(n) A.sample_fun_(n) + B);
            else
                rv = srt.disturbances.RandomVector(@(n) B.sample_fun_(n) + A);
            end
        end

        function rv = minus(A, B)
            if isa(A, 'srt.disturbances.RandomVector')
                rv = srt.disturbances.RandomVector(@(n) A.sample_fun_(n) - B);
            else
                rv = srt.disturbances.RandomVector(@(n) B.sample_fun_(n) - A);
            end
        end

        function rv = times(A, B)
            if isa(A, 'srt.disturbances.RandomVector')
                rv = srt.disturbances.RandomVector(@(n) B * A.sample_fun_(n));
            else
                rv = srt.disturbances.RandomVector(@(n) A * B.sample_fun_(n));
            end
        end

        function rv = mtimes(A, B)
            s = obj.sample();
            if isa(A, 'srt.disturbances.RandomVector')
                error(['Right matrix multiplication not supported for ', ...
                    'srt.disturbances.RandomVector']);
            else
                if size(A, 2) ~= length(s);
                    error(['Incorrect dimensions for matrix multiplication. ', ...
                        'Check that the number of columns in the matrix ', ...
                        'matches the lenght of a sample. To perform ', ...
                        'elementwise multiplication, use ''.*''.']);
                end

                rv = srt.disturbances.RandomVector(@(n) A * B.sample_fun_(n));
            end
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