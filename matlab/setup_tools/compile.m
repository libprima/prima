function compile(solvers, mexdir, fortd, gateways, options)
%COMPILE mexifies the Fortran solvers.
% solvers: list of the solvers to mexify
% mexdir: the directory that will contain the mexified solvers
% fortd: the directory containing the source files of the Fortran solvers
% gateways: the directory containing the MEX gateways of the Fortran solvers
% options: some options

% Remarks on the working directory:
% During the compilation, it is important to work in the correct directory. Otherwise, the files can
% be linked mistakenly, leading to runtime errors such as SEGFAULT.
% 1. Each [precision, debug_flag] specifies a version of the common files; they are compiled in
% work_dir = fullfile(directory_of_common_files, pdstr(precision, debug_flag));
% 2. Each [variant, precision, debug_flag] specifies a version of the solver; it is compiled in
% work_dir = fullfile(directory_of_solver_variant, pdstr(precision, debug_flag));
% 3. When compiling the solver corresponding to [variant, precision, debug_flag], we need the
% module and object files in
% common_dir = fullfile(directory_of_common_files, pdstr(precision, debug_flag)).
% 4. All the working directories (i.e., `work_dir`) should be sanitized (i.e., removing the existing
% module and object files) before the compilation.
%
% Remarks on the compilation options -O and -g:
% -O and -g may lead to (slightly) different behaviors of the mexified code. This was observed
% on 2021-09-09 in a test of NEWUOA on the AKIVA problem of CUTEst. It was because the mexified code
% produced different results when it was supposed to evaluate COS(0.59843577329095299_DP) amid OTHER
% CALCULATIONS: with -O, the result was 0.82621783366991353; with -g, it became 0.82621783366991364.
% Bizarrely, if we write a short Fortran program to evaluate only COS(0.59843577329095299_DP),
% then the result is always 0.82621783366991364, regardless of -O or -g. No idea why.

% COMPILE starts

% Directories
cpwd = pwd();  % The current directory, which may not be the directory containing this m file.
% `modern_fortd`: the directory containing the modernized source files of the Fortran solvers
modern_fortd = fortd;
% `classical_fortd`: the directory containing the classical source files of the Fortran solvers
classical_fortd = fullfile(fortd, 'classical');
% `common`: the directory that contains some common source files shared by all the Fortran solvers
common = fullfile(fortd, 'common');

% `options.debug_only` and `options.debug` indicate whether to compile the debugging version of the
% solvers. `debug_only` prevails if both of them are present (e.g., debug_only = true, debug = false).
% `debug_only` is needed only during the development to save compilation time.
if isfield(options, 'debug_only') && islogicalscalar(options.debug_only) && options.debug_only
    debug_flags = {true};
elseif isfield(options, 'debug') && islogicalscalar(options.debug) && options.debug
    debug_flags = {true, false};
else
    debug_flags = {false};  % This is the default: only compile the non-debugging optimized version.
end
precisions = all_precisions();
variants = all_variants();

% `options.verbose` indicates whether to do the compilation in the verbose mode.
verbose = (isfield(options, 'verbose') && islogicalscalar(options.verbose) && options.verbose);
if verbose
    verbose_option = '-v';
else
    verbose_option = '-silent';
end

% Name of the file that contains the list of Fortran files. There should be such a file in each
% Fortran source code directory, and the list should indicate the dependence among the files.
filelist = 'ffiles.txt';


% Compile the common files. They are shared by all solvers. We compile them only once.

% debug.F contains debugging subroutines tailored for MEX.
copyfile(fullfile(gateways, 'debug.F'), common);
% ppf.h contains preprocessing directives. It is needed only when compiling the common files.
header_file = fullfile(common, 'ppf.h');
header_file_bak = fullfile(common, 'ppf.bak');
copyfile(header_file, header_file_bak);
% Common Fortran source files.
common_files = [list_files(common, filelist), fullfile(gateways, 'fmxapi.F'), fullfile(gateways, 'cbfun.F')];

fprintf('Compiling the common files ... ');
for idbg = 1 : length(debug_flags)
    mex_options = {verbose_option, ['-', dbgstr(debug_flags{idbg})]};
    for iprc = 1 : length(precisions)
        prepare_header(header_file, precisions{iprc}, debug_flags{idbg});
        work_dir = fullfile(common, pdstr(precisions{iprc}, debug_flags{idbg}));
        prepare_work_dir(work_dir);
        cd(work_dir);
        % We can NOT write the loop below as `mex(mex_options{:}, '-c', common_files{:});`
        % Because such a command may not respect the order of common_files{:}, which is critical here.
        for icf = 1 : length(common_files)
            if verbose
                mex(mex_options{:}, '-c', common_files{icf});
            else
                evalc('mex(mex_options{:}, ''-c'', common_files{icf})');  % Suppress the output.
            end
            % The module/object files are dumped to the current directory, namely `work_dir`.
        end
    end
end
fprintf('Done.\n');


% Compile the solvers.
for isol = 1 : length(solvers)
    solver = solvers{isol};
    fprintf('Compiling %s ... ', solver);
    gateway = fullfile(gateways, [solver, '_mex.F']);
    for ivar = 1 : length(variants)
        if strcmp(variants{ivar}, 'classical')
            soldir = fullfile(classical_fortd, solver);
        else
            soldir = fullfile(modern_fortd, solver);
        end
        for idbg = 1 : length(debug_flags)
            if strcmp(variants{ivar}, 'classical') && debug_flags{idbg}
                % The support for the classical variant is limited. No debugging version.
                continue
            end
            mex_options = {verbose_option, ['-', dbgstr(debug_flags{idbg})]};
            for iprc = 1 : length(precisions)
                work_dir = fullfile(soldir, pdstr(precisions{iprc}, debug_flags{idbg}));
                prepare_work_dir(work_dir);
                common_dir = fullfile(common, pdstr(precisions{iprc}, debug_flags{idbg}));
                copyfiles(list_mod_files(common_dir), work_dir);
                cd(work_dir);
                src_files = list_files(soldir, filelist);
                % We can NOT write the loop below as `mex(mex_options{:}, '-c', common_files{:});`
                % Because such a command may not respect the order of common_files{:}, which is critical here.
                for isf = 1 : length(src_files)
                    if verbose
                        mex(mex_options{:}, '-c', src_files{isf});
                    else
                        evalc('mex(mex_options{:}, ''-c'', src_files{isf})');  % Suppress the output.
                    end
                    % The module/object files are dumped to the current directory, namely `work_dir`.
                end
                obj_files = [list_obj_files(common_dir), list_obj_files(work_dir)];
                mexname = get_mexname(solver, precisions{iprc}, debug_flags{idbg}, variants{ivar});
                if verbose
                    mex(mex_options{:}, obj_files{:}, gateway, '-output', mexname, '-outdir', mexdir);
                else
                    evalc('mex(mex_options{:}, obj_files{:}, gateway, ''-output'', mexname, ''-outdir'', mexdir)'); % Suppress the output.
                end
            end
        end
    end
    fprintf('Done.\n');
end


% Restore header_file.
if exist(header_file_bak, 'file')
    movefile(header_file_bak, header_file);
end

cd(cpwd);  % Go back to `cpwd`.


% COMPILE ends
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function files = list_files(directory, filelist)
%LIST_FILES lists the files in `directory` according to `filelist`, which should be a plain text
% file under `directory`.

files = regexp(fileread(fullfile(directory, filelist)), '\n', 'split');
files = strtrim(files(~cellfun(@isempty, files)));
files = fullfile(directory, files);

% LIST_FILES ends
return


function prepare_header(header_file, precision, debug_flag)
%PREPARE_HEADER prepares `header_file` for the compilation according to `precision` and `debug_flag`.

switch precision
case {'s', 'single'}
    rep_str(header_file, '#define REAL_PRECISION_ 64', '#define REAL_PRECISION_ 32');
    rep_str(header_file, '#define REAL_PRECISION_ 128', '#define REAL_PRECISION_ 32');
    rep_str(header_file, '#define QP_AVAILABLE_ 1', '#define QP_AVAILABLE_ 0');
case {'q', 'quadruple'}
    rep_str(header_file, '#define REAL_PRECISION_ 32', '#define REAL_PRECISION_ 128');
    rep_str(header_file, '#define REAL_PRECISION_ 64', '#define REAL_PRECISION_ 128');
    rep_str(header_file, '#define QP_AVAILABLE_ 0', '#define QP_AVAILABLE_ 1');
otherwise
    rep_str(header_file, '#define REAL_PRECISION_ 32', '#define REAL_PRECISION_ 64');
    rep_str(header_file, '#define REAL_PRECISION_ 128', '#define REAL_PRECISION_ 64');
    rep_str(header_file, '#define QP_AVAILABLE_ 1', '#define QP_AVAILABLE_ 0');
end

if debug_flag
    rep_str(header_file, '#define DEBUGGING_ 0', '#define DEBUGGING_ 1');
else
    rep_str(header_file, '#define DEBUGGING_ 1', '#define DEBUGGING_ 0');
end

% PREPARE_HEADER ends
return


function s = pdstr(precision, debug_flag)
%PDSTR returns a string according to `precision` and `debug_flag`.
s = [precision(1), dbgstr(debug_flag)];
% PDSTR ends
return


function prepare_work_dir(directory)
%PREPARE_WORKDIR prepares `directory` for the compilation: if it does not exist, create it;
% otherwise, clean it up.
if exist(directory, 'dir')
    % Clean up `directory` so that it is proper for the compilation. Without doing this, files may
    % be linked mistakenly, leading to runtime errors such as SEGFAULT.
    cellfun(@(filename) delete(filename), list_modo_files(directory));
else
    mkdir(directory);
end
% PREPARE_WORKDIR ends
return


function mod_files = list_mod_files(dir_name)
%LIST_MOD_FILES lists all module files (*.mod) in a directory

mod_files = files_with_wildcard(dir_name, '*.mod');

return


function obj_files = list_obj_files(dir_name)
%LIST_OBJ_FILES lists all object files (*.o, *.obj) in a directory

obj_files = [files_with_wildcard(dir_name, '*.o'), files_with_wildcard(dir_name, '*.obj')];

return


function modo_files = list_modo_files(dir_name)
%LIST_MODO_FILES lists all module or object files in a directory

modo_files = [list_mod_files(dir_name), list_obj_files(dir_name)];

return
