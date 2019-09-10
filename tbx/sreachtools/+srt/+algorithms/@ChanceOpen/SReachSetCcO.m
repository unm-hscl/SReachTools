function varargout = SReachSetCcO(method_str, sys, prob_thresh, safety_tube, ...
    options)


    myeps = 1e-10; % Proxy for 0. Ideally, must be eps but MPT needs slack
    validatestring(method_str,{'chance-open'});

    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('prob_thresh', @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','>=',0.5,'<=',1}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'Tube'}, ...
        {'nonempty'}));

    try
        inpar.parse(sys, prob_thresh, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % Ensure that options are provided are appropriate
    otherInputHandling(method_str, sys, options);


end

function extra_info = create_dummy_extra_info_wmax(xmax_soln)
    % Create a struct for extra_info_wmax since not computed with dummy
    % values for the opt_XXXXX and vertices_underapprox_polytope
    extra_info.xmax = xmax_soln.xmax;
    extra_info.Umax = xmax_soln.Umax;
    extra_info.xmax_reach_prob = xmax_soln.reach_prob;
    extra_info.opt_theta_i = [];
    extra_info.opt_input_vec_at_vertices = [];
    extra_info.opt_reach_prob_i = [];
    extra_info.vertices_underapprox_polytope = [];
end

function extra_info = create_dummy_extra_info_empty()
    % Create a struct for extra_info when not computed
    extra_info.xmax = [];
    extra_info.Umax = [];
    extra_info.xmax_reach_prob = [];
    extra_info.opt_theta_i = [];
    extra_info.opt_input_vec_at_vertices = [];
    extra_info.opt_reach_prob_i = [];
    extra_info.vertices_underapprox_polytope = [];
end

function otherInputHandling(method_str, sys, options)
    % Input handling for SReachSetCcO

    % Consider updating SReachSetGpO.m if any changes are made here

    % Ensure Gaussian-perturbed system
    validateattributes(sys.dist, {'RandomVector'}, {'nonempty'},...
        'SReachSetCcO/otherInputHandling', 'sys.dist');
    validatestring(sys.dist.type, {'Gaussian'}, {'nonempty'},...;
        'SReachSetCcO/otherInputHandling', 'sys.dist.type');

    % Check if prob_str and method_str are consistent
    if ~strcmpi(options.prob_str,'term')
        throwAsCaller(SrtInvalidArgsError(['Mismatch in prob_str in the ', ...
            'options']));
    end
    if ~strcmpi(options.method_str,method_str)
        throwAsCaller(SrtInvalidArgsError(['Mismatch in method_str in the ', ...
            'options']));
    end

    % Make sure the user specified set_of_dir_vecs is valid
    if size(options.set_of_dir_vecs, 1) ~= sys.state_dim
        throwAsCaller(SrtInvalidArgsError(sprintf(['set_of_dir_vecs should', ...
            ' be a collection of %d-dimensional column vectors.'], ...
            sys.state_dim)));
    end
    if size(options.set_of_dir_vecs, 2) < 2, ...
        throwAsCaller(SrtInvalidArgsError(['set_of_dir_vecs should at ', ...
            'least have two directions.']));
    end
    % Make sure the user specified init_safe_set_affine is of the correct
    % dimension
    if options.init_safe_set_affine.Dim ~= sys.state_dim &&...
            ~isempty(options.init_safe_set_affine.H)
        throwAsCaller(SrtInvalidArgsError(['init_safe_set_affine must be ', ...
            'an sys.state_dim-dimensional affine set']));
    end
end
