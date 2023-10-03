function recursive(solver, options)
%RECURSIVE verifies that the solvers can be called recursively.

if nargin < 2
    options = struct();
end

% Set the dimension of the problem
if isfield(options, 'n')
    n = options.n;
else
    if ismember(solver, {'uobyqa', 'cobyla'})
        n = 20;
    else
        n = 50;
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
    opt = struct();
    opt.verbose = false;
    opt.debug = true;
    opt.debug_only = true;
    opt.classical = false;
    opt.single = false;
    setup(solver, opt);
    cd(old_directory);
end
solver_name = solver;
solver = str2func(solver);

% Define the objective function, which is based on the Rosenbrock function and recursive calls of
% the solver.
fun = @chrosen;
for i = 1 : depth
    fun = @(x) rfun(x, fun, solver, n);
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

% Turn off the warning about the debug mode.
orig_warning_state = warning;
warning('off', [solver_name, ':Debug']);

% Conduct the test
tic;
fprintf('\n>>>>>> Recursive test for %s starts <<<<<<\n', solver_name);

% Call the solver
opt = struct();
opt.iprint = 3;
opt.debug = true;
opt.rhoend = 1.0e-3;
opt.maxfun = min(100*n, 5e3);

% We call the solver two times, in case something does not finish correctly during the first run.
[x, fx, exitflag, output] = solver(fun, randn(n, 1), opt)
[x, fx, exitflag, output] = solver(fun, randn(n, 1), opt)

fprintf('\n>>>>>> Recursive test for %s ends <<<<<<\n', solver_name);
toc;

% Restore the random number generator state
rng(orig_rng_state);
% Restore the warning state
warning(orig_warning_state);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = rfun(x, fun, solver, n)
%RFUN defines a function of x by minimizing fun([x; y]) with respect to y in R^n using a solver.
opt.debug = true;
opt.rhoend = 1.0e-3;
[~, f] = solver(@(y) fun([x; y]), randn(2, 1), opt);
return
