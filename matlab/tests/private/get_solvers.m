function spaths = get_solvers(solvers, test_dir, compile_flag)

olddir = pwd();  % Record the current directory
oldpath = path();  % Record the current path.

% We allow `solvers` to be the name of a particular solver.
if isa(solvers, 'char') || isa(solvers, 'string')
    solvers = {solvers};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ready_solvers = {'newuoa', 'cobyla'};  % Solvers whose development is (almost) finished.
test_ready_solvers = ~isempty(intersect(lower(solvers), ready_solvers));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directories.
neupdfo_dir = test_dir;  % `test_dir` is the root directly of a copy of the package made for the test.
callstack = dbstack;
invoker = callstack(2).name;  % The function that calls this function.
isverify = strcmp(invoker, 'verify');  % Are we conduction verification?
isprofile = strcmp(invoker, 'profile');  % Are we profiling the solvers?

% Get the solver paths.
pdfo_dir = fullfile(neupdfo_dir, 'OPDFO');
solver_dir = fullfile(pdfo_dir, 'matlab', 'interfaces');
solvern_dir = fullfile(neupdfo_dir, 'matlab', 'interfaces');
spaths={solver_dir, solvern_dir};

% Set up the solvers, taking particularly `compile_flag` (true/false) into account.

try
    for is = 1 : length(solvers)

        if strcmp(solvers{is}(end), 'n')
            tested_solver_name = solvers{is}(1:end-1);
        else
            tested_solver_name = solvers{is};
        end

        if compile_flag  % Compilation requested.

            % Define the compilation options.
            mexopt = struct();

            % During verification, save compilation time for non-ready solvers by specifying
            % debug_only = true: only the debugging version will be compiled.
            % N.B.: In this case, the classical variant will NOT be compiled and hence not be tested
            % even if mexopt.classical is set to true.
            mexopt.debug_only = isverify && ~test_ready_solvers;

            % When we are not in verification, only the non-debugging version will be compiled.
            mexopt.debug = isverify;

            % Save compilation time during the verification of non-ready solvers by not compiling
            % the classical variant; focus on the verification of the modernized variant.
            mexopt.classical = (isverify && test_ready_solvers) || isprofile;

            % Save compilation time during the verification of non-ready solvers by not compiling
            % the quadruple precision; focus on the verification of the double/single precision.
            mexopt.quadruple = isverify && test_ready_solvers;

            % Include the single precision into the verification.
            mexopt.single = isverify;

            cd(neupdfo_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup([tested_solver_name, 'n'], mexopt);

            if isverify
                cd(pdfo_dir);
                clear('setup');  % Without this, the next line may not call the latest version of `setup`
                setup(tested_solver_name, mexopt);
            end

        else  % No compilation. Set up the path only.

            cd(neupdfo_dir);
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
    cd(olddir);  % Go back to olddir.
    setpath(oldpath);  % Restore the path to oldpath.
    rethrow(exception);
end

cd(olddir);
