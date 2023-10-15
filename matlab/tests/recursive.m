function recursive(solver, options)
%RECURSIVE verifies that the solvers can be called recursively.

% Turn off the warning about the debug mode.
orig_warning_state = warning;
warning('off', [solver, ':Debug']);

if nargin < 2
    options = struct();
end

% Set the random seed. We ALTER THE SEED WEEKLY to test the solvers as much as possible.
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
orig_rng_state = rng();  % Save the current random number generator settings
rng(random_seed);

% Set the dimension of the problem
if isfield(options, 'n')
    n = options.n;
else
    if ismember(solver, {'uobyqa'})
        n = 5;
    else
        n = 10;
    end
end

% Set the recursion depth
if isfield(options, 'depth')
    depth = options.depth;
else
    depth = 3;
end

% Set up the solver
if ~isfield(options, 'compile') || options.compile
    old_directory = pwd();
    cd(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    compile_options = struct();
    compile_options.verbose = false;
    compile_options.debug = (rand() < 0.5);
    compile_options.classical = false;
    compile_options.single = false;
    setup(solver, compile_options);
    cd(old_directory);
end
solver_name = solver;
solver = str2func(solver);

% Define the objective function, which is based on the Rosenbrock function and recursive calls of
% the solver.
solver_options = struct();
solver_options.debug = (rand() < 0.5);
solver_options.rhoend = 1.0e-3;
solver_options.maxfun = min(100*n, 5e3);
fun = @chrosen;
for i = 1 : depth
    fun = @(x) rfun(x, fun, solver, n, solver_options);
end

% Conduct the test
tic;
fprintf('\n>>>>>> Recursive test for %s starts <<<<<<\n', solver_name);

% Call the solver
% We call the solver two times, in case something does not finish correctly during the first run.
solver_options.iprint = 3;
[x, fx, exitflag, output] = solver(fun, randn(n, 1), solver_options)
[x, fx, exitflag, output] = solver(fun, randn(n, 1), solver_options)

fprintf('\n>>>>>> Recursive test for %s ends <<<<<<\n', solver_name);
toc;

% Restore the random number generator state
rng(orig_rng_state);
% Restore the warning state
warning(orig_warning_state);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = rfun(x, fun, solver, n, solver_options)
%RFUN defines a function of x by minimizing fun([x; y]) with respect to y in R^2 using a solver.
solver_options.iprint = 0;
solver_options.rhoend = 1.0e-2;
[~, f] = solver(@(y) fun([x; y]), randn(2, 1), solver_options);
return
