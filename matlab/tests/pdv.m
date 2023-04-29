function pdv()
% PDV tests the Fortran solvers for all precision, debugging flags, and variants.

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
options.competitor = root_dir_name;
options.compile = true;
test_dir = prepare_test_dir(fake_solver_name, funname, options);

exception = [];

try

    % Go to the test directory.
    solver_dir = fullfile(test_dir, root_dir_name);
    cd(solver_dir);

    % Compile the solvers.
    clear('setup');
    opt=struct();
    opt.verbose = true;
    opt.single=true;
    opt.quadruple=true;
    opt.debug=true;
    opt.classical=true;
    tic
    setup(opt);
    testprima
    toc

    solvers = {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'};
    precisions = {'double', 'single', 'quadruple'};
    debug_flags = {true, false};
    variants = {'modern', 'classical'};

    % Show current path information.
    showpath(solvers);

    % Test the solvers.
    fun = @sin;
    x0 = 1;
    for isol = 1 : length(solvers)
        solver = str2func(solvers{isol});
        solver
        options = struct();
        for iprc = 1 : length(precisions)
            options.precision = precisions{iprc};
            for idbg = 1 : length(debug_flags)
                options.debug = debug_flags{idbg};
                for ivar = 1 : length(variants)
                    options.classical = strcmp(variants{ivar}, 'classical');
                    options.output_xhist = true;
                    options
                    format long
                    [x, f, exitflag, output] = solver(fun, x0, options)
                    if (ismember(solvers{isol}, matlab_implemented))
                        options_mat = options;
                        options_mat.fortran = false;
                        [x, f, exitflag, output] = solver(fun, x0, options_mat)
                    end
                end
            end
        end
    end

    % Show current path information again at the end of test.
    showpath(solvers);

catch exception

    % Do nothing for the moment.

end

setpath(oldpath);  % Restore the path to oldpath.
cd(olddir);  % Go back to olddir.
fprintf('\nCurrently in %s\n\n', pwd());

if ~isempty(exception)  % Rethrow any exception caught above.
    rethrow(exception);
end

end
