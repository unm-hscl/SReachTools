function [concat_tube_A, concat_tube_b] = concatenate(obj,varargin)
    %% Construction of the concatenated tube
    tube_A_mats = cell(1, length(obj));
    [tube_A_mats{:}] = obj.tube(:).A;

    tube_b_vecs = cell(1, length(obj));
    [tube_b_vecs{:}] = obj.tube(:).b;

    %% Do we send out everything or was a slice requested?
    if nargin == 1
        % Send out everything
        concat_tube_A = blkdiag(tube_A_mats{:});
        concat_tube_b = vertcat(tube_b_vecs{:});
    else
        % Send out a slice
        time_limits = varargin{1};
        validateattributes(time_limits,{'numeric'},{'vector'});
        if length(time_limits) == 2 && min(time_limits) >= 1 &&...
            max(time_limits) <= length(obj)
            if time_limits(2) >= time_limits(1)
                % a <= b => send out the appropriate slices
                concat_tube_A =...
                    blkdiag(tube_A_mats{time_limits(1):time_limits(2)});
                concat_tube_b =...
                    vertcat(tube_b_vecs{time_limits(1):time_limits(2)});
            else
                % a > b => send out empty sets
                concat_tube_A = [];
                concat_tube_b = [];
            end
        else
            % time_limits = [a,b] is not a 2x1/1x2 vector OR it does not
            % satisfy 1<= a,b <= length(obj)
            throw(SrtInvalidArgsError('Invalid time range'));
        end
    end
end
