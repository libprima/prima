function spaths = get_solvers(solvers, test_dir, options)
%GET_SOLVERS set up the solvers for testing.
% Possible members of `solvers`:
% SOLVER, a member of {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'}.
% SOLVER_norma, the version SOLVER in the ".development/norma" directory.
% SOLVER_classical|_single|_quadruple, the classical/single-precision/quadruple-precision version of SOLVER.
% SOLVER_archiva, the version of SOLVER in the "norma" directory under `dev_arch`, which is equivalent
% to the latest archiva version of SOLVER.

% We allow `solvers` to be the name of a particular solver.
if isa(solvers, 'char') || isa(solvers, 'string')
    solvers = {solvers};
end

% We do not allow `solvers` to contain both XXXX_norma and YYYY_archiva even if XXXX and YYYY are different.
assert(~(any(endsWith(solvers, '_norma')) && any(endsWith(solvers, '_archiva'))));

% `invoker` is the function that calls this function.
callstack = dbstack;
invoker = callstack(2).name;
isverify = strcmp(invoker, 'verify');  % Are we conducting verification?
isprofile = strcmp(invoker, 'profile');  % Are we profiling?

% Record the current path.
oldpath = path();
% Record the current directory
olddir = pwd();

% Directories.
% `test_dir` is the root directory of a copy of the package made for the test.
prima_dir = test_dir;
% Path for SOLVER
solver_dir = fullfile(prima_dir, 'matlab', 'interfaces');
% Path for SOLVER_norma
norma_dir = fullfile(prima_dir, '.development', 'norma');
solvern_dir = fullfile(norma_dir, 'matlab', 'interfaces');
% The following lines get the path for SOLVER_archiva
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% `dev_arch` a subdirectory of fullfile(s_root_dir, 'archiva'). It contains the "archiva" version of solvers
% used as a benchmark for the development of the current version of the solvers.
archiva_dir_name = 'dev_arch';
% Define `archiva_dir` as the "norma" directory under the `dev_arch` directory. Indeed, the solvers in
% fullfile(prima_dir, '.development', 'archiva', archiva_dir_name) and
% fullfile(prima_dir, '.development', 'archiva', archiva_dir_name, 'norma')
% are equivalent. We use the latter because the name of the solver there is SOLVER_norma, which is
% convenient for the test.
archiva_dir = fullfile(prima_dir, '.development', 'archiva', archiva_dir_name, 'norma');
solvera_dir = fullfile(archiva_dir, 'matlab', 'interfaces');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define `spaths`, a cell array of solver paths.
spaths = {solver_dir};
if any(endsWith(solvers, '_norma'))
    spaths = [spaths, {solvern_dir}];
end
if any(endsWith(solvers, '_archiva'))
    spaths = [spaths, {solvera_dir}];
end

% Will we compile the solvers?
compile_flag = ~isfield(options, 'compile') || options.compile;

% Define `mexopts`, a cell array of structures, each of which contains the options for mexifying the
% corresponding solver.
mexopts = cell(length(solvers), 1);
for is = 1 : length(solvers)
    mexopts{is} = struct();
    if compile_flag
        % Do we compile the debugged version?
        % Yes if we are in verification or if options.debug is true.
        mexopt.debug = isverify || (isfield(options, 'debug') && options.debug);

        % Do we compile ONLY the debugging version?
        % Yes if we are profiling and options.debug is true.
        mexopt.debug_only = isprofile && (isfield(options, 'debug') && options.debug);

        % Do we compile the classical version?
        % Yes if we are in verification (unless options.no_classical = true) or if the solver name
        % ends with '_classical' or SOLVER_classical is requested.
        mexopts{is}.classical = (isverify && ~(isfield(options, 'no_classical') && options.no_classical)) ...
            || endsWith(solvers{is}, '_classical') || ismember([solvers{is}, '_classical'], solvers);


        % Do we compile the single-precision version?
        % Yes if we are in verification or if the solver name ends with '_single' or SOLVER_single is requested.
        mexopts{is}.single = (isverify || endsWith(solvers{is}, '_single') || ismember([solvers{is}, '_single'], solvers));

        % Do we compile the quadruple-precision version?
        % Yes if we are in verification or if the solver name ends with '_quadruple' or SOLVER_quadruple is requested.
        mexopts{is}.quadruple = (isverify || endsWith(solvers{is}, '_quadruple') || ismember([solvers{is}, '_quadruple'], solvers));

        % Should we be verbose?
        mexopts{is}.verbose = (isfield(options, 'verbose') && options.verbose);
    end
end

% SOLVER_classical is obtained by preparing SOLVER. Thus we remove solvers ending with '_classical'.
% The same for _single and _quadruple.
[solvers, ind] = unique(regexprep(solvers, '(_classical|_single|_quadruple)', ''));
mexopts = mexopts(ind);

% Compile the solvers.
exception = [];
try

    compiler_options_modified = compile_flag && isfield(options, 'compiler_options') && ...
        (isa(options.compiler_options, 'char') || isa(options.compiler_options, 'string'));
    if compiler_options_modified
        set_compiler_options(options.compiler_options);
    end

    for is = 1 : length(solvers)

        solver = solvers{is};

        % The following `cd` decides which version of the solver to compile.
        if endsWith(solver, '_norma')
            cd(norma_dir);  % Compile SOLVER_norma
        elseif endsWith(solver, '_archiva')
            solver = regexprep(solver, '_archiva', '_norma');  % Compile SOLVER_archiva
            cd(archiva_dir);
        else  % SOLVER or SOLVER_classical|_single|_quadruple
            cd(prima_dir);  % Compile SOLVER
        end

        clear('setup');  % Without this, the next line may not call the latest version of `setup`

        if compile_flag  % Compilation requested.
            setup(solver, mexopts{is});
        else  % No compilation. Set up the path only.
            setup('path');
        end

    end

catch exception
    % Do nothing for the moment.
end

% Restore the compiler options.
if compiler_options_modified
    restore_compiler_options();
end
% Go back to the old directory.
cd(olddir);

% If there is an exception, restore the path and rethrow the exception.
if ~isempty(exception)
    setpath(oldpath);
    rethrow(exception);
end
