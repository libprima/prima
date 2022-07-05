function pdv
% PDV tests the Fortran solvers for all precision, debugging flags, and variants.

oldpath = path();  % Record the current path.
restoredefaultpath;  % Restore the "right out of the box" path of MATLAB

olddir = pwd();  % Record the current directory.

% Prepare the test directory, i.e., `test_dir`.
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
fake_solver_name = funname;
options.compile = true;
test_dir = prepare_test_dir(fake_solver_name, funname, options);

exception = [];

try

    % Go to the test directory. This is not really necessary. It will not affect the test, but any
    % output (e.g., NEWUOA_output.txt, fort.6) will be dumped to `test_dir`.
    cd(test_dir);

    % Compile the solvers.
    clear('setup');
    opt=struct();
    opt.single=true;
    opt.quadruple=true;
    opt.debug=true;
    opt.classical=true;
    tic
    setup(opt);
    testpdfon
    toc

    solvers = {'cobylan', 'uobyqan', 'newuoan', 'bobyqan', 'lincoan'};
    %solvers = {'cobylan', 'newuoan', 'bobyqan', 'lincoan'};  % uobyqan cannot solve 1D problem as of 20220502
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
