function pdv
% PDV tests the Fortran solvers for all precision, debugging flags, and variants.

oldpath = path();  % Record the current path.
restoredefaultpath;  % Restore the "right out of the box" path of MATLAB

olddir = pwd();  % Record the current directory.

% Prepare the test directory, i.e., `test_dir`.
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
test_dir = prepare_test_dir(funname);

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
    testpdfo
    toc

    % Show current path information.
    showpath();

    % Test the solvers.
    solvers = {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'};
    precisions = {'double', 'single', 'quadruple'};
    debug_flags = {true, false};
    variants = {'modern', 'classical'};
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
                    options
                    format long
                    [a, b, c, d] = solver(fun, x0, options)
                end
            end
        end
    end

catch exception

    % Do nothing for the moment.

end

setpath(oldpath);  % Restore the path to oldpath.
cd(olddir);  % Go back to olddir.

if ~isempty(exception)  % Rethrow any exception caught above.
    rethrow(exception);
end

end
