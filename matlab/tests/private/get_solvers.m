function spaths = get_solvers(solvers, test_dir, options)
%GET_SOLVERS set up the solvers for testing. `spaths` will be a cell array containing paths to
% these solvers.
% Possible members of `solvers`:
% SOLVER, a member of {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'}.
% SOLVER_norma, the version SOLVER in the ".development/norma" directory.
% SOLVER_classical|_default|_half|_single|_quadruple, the classical/default/single-precision/quadruple-precision version of SOLVER.
% SOLVER_archiva, the version of SOLVER in the "norma" directory under `dev_arch`, which is equivalent
% to the latest archiva version of SOLVER.

spaths = cell(1, length(solvers));

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

% mfilepath: directory where this .m file resides.
mfilepath = fileparts(mfilename('fullpath'));

% root_dir: root directory of the project
root_dir = fileparts(fileparts(fileparts(mfilepath)));
[~, root_dir_name] = fileparts(root_dir);

% Record the current path.
oldpath = path();
% Record the current directory
olddir = pwd();

% Will we compile the solvers?
compile_flag = ~isfield(options, 'compile') || options.compile;

% Do we have compiler options?
with_compiler_options = compile_flag && isfield(options, 'compiler_options') && ischarstr(options.compiler_options);

% Define `mexopts`, a cell array of structures, each of which contains the options for mexifying the
% corresponding solver.
mexopts = cell(length(solvers), 1);
for is = 1 : length(solvers)
    mexopts{is} = struct();
    if compile_flag
        % Do we compile the debug version?
        % Yes if we are in verification or if options.debug is true.
        mexopts{is}.debug = isverify || (isfield(options, 'debug') && options.debug);

        % Do we compile ONLY the debugging version?
        % Yes if we are profiling and options.debug is true.
        %mexopts{is}.debug_only = isprofile && (isfield(options, 'debug') && options.debug);

        % Do we compile the classical version?
        % Yes if we are in verification (unless options.no_classical = true) or if the solver name
        % ends with '_classical' or SOLVER_classical is requested.
        mexopts{is}.classical = (isverify && ~(isfield(options, 'no_classical') && options.no_classical)) ...
            || endsWith(solvers{is}, '_classical') || ismember([solvers{is}, '_classical'], solvers);

        % Do we compile the half-precision version?
        % Yes if we are in verification or if the solver name ends with '_half' or
        % SOLVER_half is requested or options.precision='half'.
        mexopts{is}.half = (isverify || endsWith(solvers{is}, '_half') || ismember([solvers{is}, '_half'], solvers)) || ...
            (isfield(options, 'precision') && (isa(options.precision, 'char') ||  isa(options.precision, 'string')) && strcmpi(options.precision, 'half'));
        mexopts{is}.half = mexopts{is}.half && half_precision_available();

        % Do we compile the single-precision version?
        % Yes if we are in verification or if the solver name ends with '_single' or
        % SOLVER_single is requested or options.precision='single'.
        mexopts{is}.single = (isverify || endsWith(solvers{is}, '_single') || ismember([solvers{is}, '_single'], solvers)) || ...
            (isfield(options, 'precision') && (isa(options.precision, 'char') ||  isa(options.precision, 'string')) && strcmpi(options.precision, 'single'));

        % Do we compile the quadruple-precision version?
        % Yes if we are in verification or if the solver name ends with '_quadruple' or
        % SOLVER_quadruple is requested or options.precision is 'quadruple'.
        mexopts{is}.quadruple = (isverify || endsWith(solvers{is}, '_quadruple') || ismember([solvers{is}, '_quadruple'], solvers)) || ...
            (isfield(options, 'precision') && (isa(options.precision, 'char') ||  isa(options.precision, 'string')) && strcmpi(options.precision, 'quadruple'));


        % Should we be verbose?
        mexopts{is}.verbose = (isfield(options, 'verbose') && options.verbose || ...
            (with_compiler_options && ~(isfield(options, 'verbose') && ~options.verbose && ~isprofile)));
    end
end

% SOLVER_classical is obtained by preparing SOLVER. Thus we remove solvers ending with '_classical'.
% The same for _default, _half, _single, and _quadruple.
[solvers, ind] = unique(regexprep(solvers, '(_classical|_default|_half|_single|_quadruple)', ''));
mexopts = mexopts(ind);

% Compile the solvers.
exception = [];
try

    if with_compiler_options
        set_compiler_options(options.compiler_options);
    end

    for is = 1 : length(solvers)

        solver = solvers{is};
        fprintf('\nSetting up %s ...\n', solver);

        if endsWith(solver, '_norma')
            solver_dir = fullfile(test_dir, 'norma');
        elseif endsWith(solver, '_archiva')
            solver_dir = fullfile(test_dir, 'archiva');
            % The archiva solver name is SOLVER_norma. See the comments on archiva_dir in
            % prepare_test_dir.m for details.
            solver = regexprep(solver, '_archiva', '_norma');
        else  % SOLVER or SOLVER_classical|_default|_half|_single|_quadruple
            solver_dir = fullfile(test_dir, root_dir_name);
        end

        spaths{is} = fullfile(solver_dir, 'matlab', 'interfaces');

        clear('setup');  % Without this, the next line may not call the latest version of `setup`

        cd(solver_dir);
        fprintf('\nDirectory changed to: %s\n', pwd());


        if compile_flag  % Compilation requested.

            % Revise the header file to use the intended integer kind.
            header_file = fullfile(cd(), 'fortran', 'common', 'ppf.h');
            integer_kind = 0;
            if isfield(options, 'integer_kind')
                integer_kind = options.integer_kind;
            end
            if integer_kind ~= 0
                % Do not change the integer kind if it is 0!!! This enables us to test two solvers
                % with different integer kinds by modifying their respective header files manually.
                rep_str(header_file, '#define PRIMA_INTEGER_KIND 0', ['#define PRIMA_INTEGER_KIND ', num2str(integer_kind)]);
                rep_str(header_file, '#define PRIMA_INTEGER_KIND 16', ['#define PRIMA_INTEGER_KIND ', num2str(integer_kind)]);
                rep_str(header_file, '#define PRIMA_INTEGER_KIND 32', ['#define PRIMA_INTEGER_KIND ', num2str(integer_kind)]);
                rep_str(header_file, '#define PRIMA_INTEGER_KIND 64', ['#define PRIMA_INTEGER_KIND ', num2str(integer_kind)]);
            end

            setup(solver, mexopts{is});
        else  % No compilation. Set up the path only.
            setup('path');
        end

    end

catch exception
    % Do nothing for the moment.
end

% Restore the compiler options.
if with_compiler_options
    restore_compiler_options();
end
% Go back to the old directory.
cd(olddir);

% If there is an exception, restore the path and rethrow the exception.
if ~isempty(exception)
    setpath(oldpath);
    rethrow(exception);
end
