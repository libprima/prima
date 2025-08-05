function pdv(options)
% PDV tests the Fortran solvers for all Precision, Debugging flags, and Variants.

if nargin < 1
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
fprintf('\nThe seed is\t\t%d\n\n', yw);
% Define the random seed by yw
random_seed = yw;
orig_rng_state = rng();  % Save the current random number generator settings
rng(random_seed);

oldpath = path();  % Record the current path.
restoredefaultpath;  % Restore the "right out of the box" path of MATLAB

olddir = pwd();  % Record the current directory.

matlab_implemented = {'newuoa'};  % Solvers that has a MATLAB implementation.

% mfilepath: directory where this file resides
mfilepath = fileparts(mfilename('fullpath'));

% root_dir: root directory of the project
root_dir = fileparts(fileparts(mfilepath));
[~, root_dir_name] = fileparts(root_dir);

% Prepare the test directory, i.e., `test_dir`.
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
fake_solver_name = root_dir_name;
opt.competitor = root_dir_name;
opt.compile = true;
test_dir = prepare_test_dir(fake_solver_name, funname, opt);

exception = [];

try

    % Go to the test directory.
    solver_dir = fullfile(test_dir, root_dir_name);
    cd(solver_dir);

    % Compile the solvers.
    clear('setup');
    opt=struct();
    opt.verbose = true;
    opt.half = half_precision_available();
    opt.single = true;
    opt.double = true;
    opt.quadruple = true;
    opt.debug = true;
    %opt.classical = true;
    opt.classical = false;

    tic

    setup(opt);

    solvers = {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'};
    precisions = {'half', 'single', 'double', 'quadruple'};
    precisions = precisions([opt.half, opt.single, opt.double, opt.quadruple]);
    debug_flags = {true, false};
    %variants = {'modern', 'classical'};
    variants = {'modern'};

    % Show current path information.
    showpath(solvers);

    % Test the solvers.
    %fun = @sin; x0 = -1;
    fun = @chrosen; x0 = [-1, -1];
    x0 = x0 + 0.1*randn(size(x0));
    for isol = 1 : length(solvers)
        solver = str2func(solvers{isol});
        solver
        opt = struct();
        for iprc = 1 : length(precisions)
            opt.precision = precisions{iprc};
            for idbg = 1 : length(debug_flags)
                opt.debug = debug_flags{idbg};
                for ivar = 1 : length(variants)
                    opt.classical = strcmp(variants{ivar}, 'classical');
                    if ismac && strcmp(func2str(solver), 'cobyla') && strcmp(opt.precision, 'half') && opt.classical
                        % Skip the classical cobyla in half precision on macOS, as it will encounter an infinite cycling.
                        continue;
                    end
                    if ismac && strcmp(func2str(solver), 'bobyqa') && strcmp(opt.precision, 'quadruple') && opt.classical
                        % Skip the classical bobyqa in quadruple precision on macOS, as it will encounter a segmentation fault.
                        continue;
                    end
                    opt.output_xhist = true;
                    opt.maxfun = 200*length(x0);
                    opt.rhoend = 1.0e-4;
                    opt.iprint = randi([-4, 4]);
                    format long
                    opt
                    [x, f, exitflag, output] = solver(fun, x0, opt)
                    if (ismember(solvers{isol}, matlab_implemented))
                        opt_mat = opt;
                        opt_mat.fortran = false;
                        [x, f, exitflag, output] = solver(fun, x0, opt_mat)
                    end
                end
            end
        end
    end

    toc

    % Show current path information again at the end of test.
    showpath(solvers);

catch exception

    % Do nothing for the moment.

end

% Restore the random number generator state
rng(orig_rng_state);
setpath(oldpath);  % Restore the path to oldpath.
cd(olddir);  % Go back to olddir.
fprintf('\nCurrently in %s\n\n', pwd());

if ~isempty(exception)  % Rethrow any exception caught above.
    rethrow(exception);
end

end
