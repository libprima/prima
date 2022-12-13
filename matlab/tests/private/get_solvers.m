function spaths = get_solvers(solvers, test_dir, options)

olddir = pwd();  % Record the current directory
oldpath = path();  % Record the current path.

% We allow `solvers` to be the name of a particular solver.
if isa(solvers, 'char') || isa(solvers, 'string')
    solvers = {solvers};
end

% Directories.
prima_dir = test_dir;  % `test_dir` is the root directory of a copy of the package made for the test.
callstack = dbstack;
invoker = callstack(2).name;  % The function that calls this function.
isverify = strcmp(invoker, 'verify');  % Are we conduction verification?
isprofile = strcmp(invoker, 'profile');  % Are we profiling the solvers?

% Get the solver paths.
pdfo_dir = fullfile(prima_dir, 'OPDFO');
solver_dir = fullfile(pdfo_dir, 'matlab', 'interfaces');
solvern_dir = fullfile(prima_dir, 'matlab', 'interfaces');
spaths={solver_dir, solvern_dir};

% Set up the solvers, taking particularly `compile_flag` and `debug_flag` (true/false) into account.
compile_flag = ~isfield(options, 'compile') || options.compile;
debug_flag = isverify || (isfield(options, 'debug') && options.debug);

try

    compiler_options_modified = false;
    if isfield(options, 'compiler_options') && (isa(options.compiler_options, 'char') ...
            || isa(options.compiler_options, 'string'))
        set_compiler_options(options.compiler_options);
        compiler_options_modified = true;
    end

    for is = 1 : length(solvers)

        if strcmp(solvers{is}(end), 'n')
            tested_solver_name = solvers{is}(1:end-1);
        else
            tested_solver_name = solvers{is};
        end

        if compile_flag  % Compilation requested.

            % Define the compilation options.
            mexopt = struct();
            if isfield(options, 'verbose')
                mexopt.verbose = options.verbose;
            end

            % When we are not in verification, only the non-debugging version will be compiled.
            % We test both the debugging and non-debugging version during the verification.
            mexopt.debug_only = false;
            mexopt.debug = debug_flag;

            mexopt.classical = true;

            % Include the single and quadruple precisions into the verification.
            mexopt.single = isverify;
            mexopt.quadruple = isverify;

            cd(prima_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup([tested_solver_name, 'n'], mexopt);

            if isverify
                cd(pdfo_dir);
                clear('setup');  % Without this, the next line may not call the latest version of `setup`
                setup(tested_solver_name, mexopt);
            end

        else  % No compilation. Set up the path only.

            cd(prima_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup('path');

            if isverify
                cd(pdfo_dir);
                clear('setup');  % Without this, the next line may not call the latest version of `setup`
                setup('path');
            end

        end

    end
catch exception
    if compiler_options_modified
        restore_compiler_options();  % Restore the compiler options.
    end
    cd(olddir);  % Go back to olddir.
    setpath(oldpath);  % Restore the path to oldpath.
    rethrow(exception);
end

if compiler_options_modified
    restore_compiler_options();  % Restore the compiler options.
end
cd(olddir);
