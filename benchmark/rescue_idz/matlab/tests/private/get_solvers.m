function spaths = get_solvers(solvers, test_dir, options)
%GET_SOLVERS set up the solvers for testing.
% Possible members of `solvers`:
% SOLVER, a member of {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'}.
% SOLVER_last, the version SOLVER in the `last` directory.
% SOLVER_classical|_single|_quadruple, the classical/single-precision/quadruple-precision version of SOLVER.
% SOLVER_archiva, the version of SOLVER in the `archiva` directory under `last_xxxx`, which is equivalent
% to the latest archiva version of SOLVER.

% We allow `solvers` to be the name of a particular solver.
if isa(solvers, 'char') || isa(solvers, 'string')
    solvers = {solvers};
end

% We do not allow `solvers` to contain both XXXX_last and YYYY_archiva even if XXXX and YYYY are different.
assert(~(any(endsWith(solvers, '_last')) && any(endsWith(solvers, '_archiva'))));

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
% Path for SOLVER_last
last_dir = fullfile(prima_dir, 'last');
solverl_dir = fullfile(last_dir, 'matlab', 'interfaces');
% The following lines get the path for SOLVER_archiva
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
archiva_list = dir(fullfile(prima_dir, 'archiva'));  % List of all the archiva directories.
% The following line keeps only the directories with names that are numbers.
archiva_list = archiva_list([archiva_list.isdir] & ~isnan(str2double({archiva_list.name})));
[~, ind] = sort(str2double({archiva_list.name}), 'ascend');  % Sort the archiva directories by name.
archiva_dir_name = archiva_list(ind(end)).name;  % The name of the latest archiva directory.
% Define `archiva_dir` as the `last` directory under the latest archiva directory. Indeed, the solvers in
% fullfile(prima_dir, 'archiva', archiva_dir_name) and fullfile(prima_dir, 'archiva', archiva_dir_name, 'last')
% are equivalent. We use the latter because the name of the solver there is SOLVER_last, which is
% convenient for the test.
archiva_dir = fullfile(prima_dir, 'archiva', archiva_dir_name, 'last');
solvera_dir = fullfile(archiva_dir, 'matlab', 'interfaces');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define `spaths`, a cell array of solver paths.
spaths = {solver_dir};
if any(endsWith(solvers, '_last'))
    spaths = [spaths, {solverl_dir}];
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
        mexopts{is}.classical = mexopts{is}.classical || (isfield(options, 'classical') && options.classical);  % FOR TESTING !!!!!!

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
        if endsWith(solver, '_last')
            cd(last_dir);  % Compile SOLVER_last
        elseif endsWith(solver, '_archiva')
            solver = regexprep(solver, '_archiva', '_last');  % Compile SOLVER_archiva
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
