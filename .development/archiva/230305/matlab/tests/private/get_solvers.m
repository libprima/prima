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
isverify = strcmp(invoker, 'verify');  % Are we conducting verification?

% Get the solver paths.
norma_dir = fullfile(prima_dir, 'norma');
solverl_dir = fullfile(norma_dir, 'matlab', 'interfaces');
solver_dir = fullfile(prima_dir, 'matlab', 'interfaces');
spaths = {solverl_dir, solver_dir};

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

        tested_solver_name = regexprep(solvers{is}, '_norma', '');

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

            % We test the classical version unless otherwise specified.
            mexopt.classical = ~(isfield(options, 'no_classical') && options.no_classical);

            % Include the single and quadruple precisions into the verification.
            mexopt.single = isverify;
            mexopt.quadruple = isverify;

            cd(prima_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup(tested_solver_name, mexopt);

            if isverify
                cd(norma_dir);
                clear('setup');  % Without this, the next line may not call the latest version of `setup`
                setup([tested_solver_name, '_norma'], mexopt);
            end

        else  % No compilation. Set up the path only.

            cd(prima_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup('path');

            if isverify
                cd(norma_dir);
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
