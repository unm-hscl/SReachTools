function addsandbox(varargin)

    p = inputParser();

    addParameter(p, 'Verbose', 0, @(x) x == 0 || x == 1);

    parse(p, varargin{:});
    verbose = p.Results.Verbose;

    if verbose > 0
        fprintf('Adding SReachTools to path: developer mode\n');
    end

    thisFolder = fileparts(mfilename('fullpath'));

    addpath(fullfile(thisFolder, 'tbx'));
    addpath(fullfile(thisFolder, 'tbx', 'sreachtools'));
    addpath(fullfile(thisFolder, 'tbx', 'doc'));

end