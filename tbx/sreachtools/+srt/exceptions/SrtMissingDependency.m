classdef SrtMissingDependency < SrtBaseException
% Custom exception object for SReachTools setup errors
% ============================================================================
%
% Customized class for generating SReachTools setup errors, subclass of the
% standard MATLAB SrtBaseException class
%
% Usage:
% ------
% exc = SrtMissingDependency('error message')
%
% ============================================================================
%
% See also MException
%
% ============================================================================
%
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
%

    properties (Constant, Access = private)
        mnemonic = 'setupError';
    end

    methods
        function obj = SrtMissingDependency(varargin)
            obj@SrtBaseException(SrtMissingDependency.mnemonic, varargin{:});
       end
    end

    methods (Static)
        function id = getErrorId()
            id = [SrtBaseException.getErrorComponent(), ':', ...
                SrtInvalidArgsError.mnemonic];
        end
    end
end
