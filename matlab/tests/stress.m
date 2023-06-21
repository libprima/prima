function stress(solver, options)
% STRESS  Stress test for the solver on problems large dimensions.

if nargin < 2
    options = struct();
end

% Set up the solver
old_directory = pwd();
cd(fileparts(fileparts(fileparts(mfilename('fullpath')))));
opt.debug = true;
setup(solver, opt);
cd(old_directory);
solver_name = solver;
solver = str2func(solver);

% Whether to conduct a TOUGH test
tough_test = isfield(options, 'tough') && options.tough;

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
    if tough_test
        switch solver_name
        case 'uobyqa'
            n = 200; %120;
        case 'newuoa'
            n = 1500; %1000;
        case 'bobyqa'
            n = 1500; %1000;
        case 'lincoa'
            n = 400; %500;
        case 'cobyla'
            n = 300; %250; %400;
        end
    else
        switch solver_name
        case 'uobyqa'
            n = 200; %120;
        case 'newuoa'
            n = 800; %500;
        case 'bobyqa'
            n = 800; %500;
        case 'lincoa'
            n = 400; %300; %500;
        case 'cobyla'
            n = 150; %250;
        end
    end
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
test_options.maxfun = 200 * n;
test_options.rhobeg = 1;
test_options.rhoend = 1.0e-7;  % In this test, the solvers normally exit with RHO reaching RHOEND, except for UOBYQA
test_options.iprint = 2;
test_options.debug = true;

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

if ~isempty(exception)
    rethrow(exception);
end

if tough_test
    fprintf('\n>>>>>> TOUGH test ends <<<<<<\n\n');
else
    fprintf('\n>>>>>> Test ends <<<<<<\n\n');
end


return
