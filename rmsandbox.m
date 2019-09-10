function rmsandbox(varargin)

    p = inputParser();

    addParameter(p, 'Verbose', 0, @(x) x == 0 || x == 1);

    parse(p, varargin{:});
    verbose = p.Results.Verbose;

    if verbose > 0
        fprintf('Removing SReachTools from path \n');
    end

    thisFolder = fileparts(mfilename('fullpath'));

    rmpath(fullfile(thisFolder, 'tbx'));
    rmpath(fullfile(thisFolder, 'tbx', 'sreachtools'));
    rmpath(fullfile(thisFolder, 'tbx', 'doc'));

end