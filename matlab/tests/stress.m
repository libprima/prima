function stress(solver, options)
%STRESS  Stress test for the solver on problems large dimensions.
% N.B.: When the dimension beyond some limit, the MEX function will crash due to memory violations.
% As we can see below in the comments, the limit is much higher on Linux than on macOS and Windows.
% This is probably because MEX uses gfortran on Linux, and it puts automatic arrays on the heap by
% default, where as MEX uses ifort on macOS and Windows, and it puts automatic arrays on the stack
% by default.

% Turn off unwanted warnings
orig_warning_state = warnoff({solver});

if nargin < 2
    options = struct();
end

% Whether to conduct a TOUGH test
tough_test = isfield(options, 'tough') && options.tough;

% Which precision to test
if isfield(options, 'precision') && ischarstr(options.precision)
    precision = options.precision;
else
    precision = 'double';
end

% Set up the solver
old_directory = pwd();
cd(fileparts(fileparts(fileparts(mfilename('fullpath')))));
opt.verbose = true;
opt.debug = true;
opt.single = strcmpi(precision, 'single');
opt.quadruple = strcmpi(precision, 'quadruple');
setup(solver, opt);
cd(old_directory);
solver_name = solver;
solver = str2func(solver);

% Set the random seed using solver name. We ALTER THE SEED weekly to test the solvers as much as possible.
if isfield(options, 'yw')
    yw = options.yw;
elseif isfield(options, 'seed')
    yw = options.seed;
else
    yw = year_week('Asia/Shanghai');
end
fprintf('\nYW = %d\n', yw);
% Define the random seed by yw
random_seed = yw;

% Set the dimension of the problem
if isfield(options, 'n')
    n = options.n;
else
    %if ismac
    %    switch solver_name
    %    case 'uobyqa'
    %        n = 50; % UOBYQA crashes if n >= 60
    %    case 'newuoa'
    %        n = 500; % NEWUOA crashes if n >= 600
    %    case 'bobyqa'
    %        n = 500; % BOBYQA crashes if n >= 600
    %    case 'lincoa'
    %        n = 400; % LINCOA crashes if n >= 500
    %    case 'cobyla'
    %        n = 300; % COBYLA crashes if n >= 400
    %    end
    %elseif ispc
    %    switch solver_name
    %    case 'uobyqa'
    %        n = 60; % UOBYQA crashes if n >= 70
    %    case 'newuoa'
    %        n = 600; % NEWUOA crashes if n >= 700
    %    case 'bobyqa'
    %        n = 600; % BOBYQA crashes if n >= 700
    %    case 'lincoa'
    %        n = 500; % LINCOA crashes if n >= 600
    %    case 'cobyla'
    %        n = 400; % COBYLA crashes if n >= 500
    %    end
    %else
        if tough_test
            switch solver_name
            case 'uobyqa'
                n = 150; % UOBYQA crashes if n > 200
            case 'newuoa'
                n = 1600;
            case 'bobyqa'
                n = 1600;
            case 'lincoa'
                n = 1000;
            case 'cobyla'
                n = 800;
            end
        else
            switch solver_name
            case 'uobyqa'
                n = 100; % UOBYQA crashes if n > 200
            case 'newuoa'
                n = 800;
            case 'bobyqa'
                n = 800;
            case 'lincoa'
                n = 500;
            case 'cobyla'
                n = 400;
            end
        end
    %end
end

% Set the type of the problem
switch solver_name
case {'uobyqa', 'newuoa'}
    problem_type = 'u';
case 'bobyqa'
    problem_type = 'b';
case 'lincoa'
    problem_type = 'l';
case 'cobyla'
    problem_type = 'n';
end

% Set the options for the test
test_options = struct();
test_options.precision = precision;
test_options.maxfun = 500 * n;
test_options.rhobeg = 1;
test_options.rhoend = 1.0e-7;
test_options.iprint = 2;
%test_options.debug = true;
test_options.debug = ~(ispc || ismac);

% Generate the problem
problem = stress_problem(n, problem_type, random_seed);
problem.options = test_options;
original_problem = problem;
if tough_test
    problem = tough(original_problem, random_seed);
end

% Conduct the test
if tough_test
    fprintf('\n>>>>>> TOUGH test starts <<<<<<\n');
else
    fprintf('\n>>>>>> Test starts <<<<<<\n');
end

tic;
% For TOUGH tests, cobyla raises an error if the constraint evaluation fails at the starting point.
% In that case, we modify the random seed and redo the test.
redo = true;
while redo
    exception = [];
    try
        solver(problem);
    catch exception
    end
    if isempty(exception) || ~tough_test || ~strcmp(solver_name, 'cobyla')
        redo = false;
    else
        redo = strcmp(exception.identifier, 'cobyla:ConstraintFailureAtX0');
        if redo
            random_seed = random_seed + 1;
            problem = tough(original_problem, random_seed);
        end
    end
end
toc;

% Restore the behavior of displaying warnings
warning(orig_warning_state);

if ~isempty(exception)
    rethrow(exception);
end

if tough_test
    fprintf('\n>>>>>> TOUGH test ends <<<<<<\n\n');
else
    fprintf('\n>>>>>> Test ends <<<<<<\n\n');
end


return
