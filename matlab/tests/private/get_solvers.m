function spaths = get_solvers(solvers, test_dir, compile_flag)

olddir = pwd();  % Record the current directory
oldpath = path();  % Record the current path.

%% All possible solvers.
%known_solvers = {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'};

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
if strcmp(invoker, 'verify')
    pdfo_dir = fullfile(neupdfo_dir, 'OPDFO');
else
    pdfo_dir = fullfile(neupdfo_dir, 'PDFO');
end

% Get the solver paths.
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

            mexopt = struct();
            % Do not remove `~test_ready_solvers` in `mexopt.debug_only`; otherwise, we will be
            % comparing classical solvers and modernized ones!!!
            mexopt.debug_only = strcmp(invoker, 'verify') && ~test_ready_solvers;
            mexopt.debug = strcmp(invoker, 'verify');
            mexopt.classical = test_ready_solvers;
            mexopt.single = test_ready_solvers;
            mexopt.quadruple = test_ready_solvers;

            cd(neupdfo_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup([tested_solver_name, 'n'], mexopt);

            cd(pdfo_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup(tested_solver_name, mexopt);

        else  % No compilation. Set up the path only.

            cd(neupdfo_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup('path');

            cd(pdfo_dir);
            clear('setup');  % Without this, the next line may not call the latest version of `setup`
            setup('path');

        end

    end
catch exception
    cd(olddir);  % Go back to olddir.
    setpath(oldpath);  % Restore the path to oldpath.
    rethrow(exception);
end

cd(olddir);
