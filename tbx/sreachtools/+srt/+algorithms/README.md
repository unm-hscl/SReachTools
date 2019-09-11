# Algorithms

To implement a new algorithm, create a new class folder in this directory with a
short, but descriptive name of the algorithm.

Use the following directory structure for the algorithm:

```
+algorithms/
├── @Algorithm/
├── @ChanceAffine/
├── ...
└── @YourAlgorithm/
    ├── private/
    │   ├── validatedependencies.m
    │   ├── validateproblem.m
    │   └── validatesystem.m
    ├── YourAlgorithm.m
    ├── compute_fwd.m
    ├── compute_point.m
    └── compute_set.m
```

## Algorithm Class

Algorithms are implemented as a class, so that parameters and options can be set
on the class before the algorithm is run.

It is recommended to use an `inputParser` with name/value parameters to parse
algorithm options.

An example of a basic algorithm called `YourAlgorithm` is shown below:

```Matlab
% Algorithms must inherit from the abstract 'Algorithm' base class.
classdef YourAlgorithm < srt.algorithms.Algorithm

  properties
    prop1
  end

  methods
    function obj = YourAlgorithm(varargin)
      % YOURALGORITHM Construct an instance of the algorithm.

      % Call the parent constructor.
      % This validates the dependencies when the algorithm is declared and sets
      % certain properties such as verbosity, etc.
      obj = obj@srt.algorithms.Algorithm(varargin{:});

      p = inputParser;
      addParameter(p, '<name>', <default>, <validation_fn>);
      parse(p, varargin{:});

      obj.prop1 = p.Results.<name>;

    end
  end

end
```

Beyond these basic requirements, you can add any class methods or properties to
the algorithm class that you would like.

## Compute Functions

At runtime, the problem and system are passed to the algorithm's `compute_x.m`
method, where `x` is one of: `fwd`, `point`, or `set`.

The function signatures for the compute functions are:

```Matlab
function results = compute_fwd(obj, prb, sys, x0, varargin)
```

```Matlab
function results = compute_point(obj, prb, sys, x0, varargin)
```

```Matlab
function results = compute_set(obj, prb, sys, varargin)
```

Do not alter the function signatures when you implement your algorithm. Any
additional parameters should be set as properties of the algorithm class.

### Results

Results are usually passed back to the caller as a `struct`. There are no
requirements or limitations on the struct properties. However, it is recommended
to implement these properties in a consistent manner to avoid confusion.

## Private Validation Functions

The algorithms are set up to perform validation to ensure that dependencies
exist, and that the system and problem are compatible with the algorithm.

The following section describes how to implement dependency, problem, and system validation.

* `validatedependencies.m`
* `validateproblem.m`
* `validatesystem.m`

These functions are optional. If your algorithm has no external dependencies,
for example, it is not necessary to implement that file in your algorithm file
structure.

### Dependencies

If your algorithm has external dependencies, modify the `validatedependencies.m`
file in the `private` folder to check for the dependencies. If any of the
dependencies are unmet, the function should throw an error.

```Matlab
function validatedependencies(obj)
% VALIDATEDEPENDENCIES Checks if algorithm dependencies are present/valid.

% Requires CVX.

% Throw an error if CVX is not installed.

end
```

### Problem Validation

If your algorithm cannot handle a certain class of problems, modify the `validateproblem.m` file in the private folder.

```Matlab
function validateproblem(obj, prob)
% VALIDATEPROBLEM Validate problem.

validateattributes(problem, {'Problem'}, {'nonempty'});

if ~ismember(class(problem), {'TerminalHitting'})
  error('Problem is not supported.');
end

end
```

### System Validation

If your algorithm can only handle a certain class of stochastic system, or if the system disturbance must have a particular distribution, modify the `validatesystem.m` file in the private folder.

```Matlab
function validatesystem(obj, sys)
% VALIDATESYSTEM Validate system.

validateattributes(sys, {'LTISystem', 'LTVSystem'}, {'nonempty'});

end
```
