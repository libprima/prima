function compile(solvers, mexdir, modern_src, classical_src, common, gateways, options)
%COMPILE mexifies the Fortran solvers.
% solvers: list of the solvers to mexify
% mexdir: the directory that will contain the mexified solvers
% modern_src: the directory containing the modernized source files of the Fortran solvers
% classical_src: the directory containing the classical source files of the Fortran solvers
% common: the directory that contains some common source files shared by all the Fortran solvers
% gateways: the directory containing the MEX gateways of the Fortran solvers
% options: some options

% Remarks on the working directory:
% During the compilation, it is important to work in the correct directory. Otherwise, the files can
% be linked mistakenly, leading to runtime errors such as SEGFAULT.
% 1. Each [precision, debug_flag] specifies a version of the common files; they are compiled in
% work_dir = fullfile(common, pdstr(precision, debug_flag));
% 2. Each [variant, precision, debug_flag] specifies a version of the solver; it is compiled in
% work_dir = fullfile(directory_of_solver_variant, pdstr(precision, debug_flag));
% 3. When compiling the solver corresponding to [variant, precision, debug_flag], we need the
% module and object files in
% common_dir = fullfile(common, pdstr(precision, debug_flag)).
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

cpwd = pwd();

% `options.debug` indicates whether to compile the debugging version of the solvers.
if isfield(options, 'debug') && islogicalscalar(options.debug) && options.debug
    debug_flags = {true, false};
else
    debug_flags = {false};
end
precisions = all_precisions();
variants = all_variants();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ready_solvers = {'cobylan', 'newuoan'};  % To be removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Name of the file that contains the list of Fortran files. There should be such a file in each
% Fortran source code directory, and the list should indicate the dependence among the files.
filelist = 'ffiles.txt';

% Detect whether we are running a 32-bit MATLAB, where maxArrayDim = 2^31-1, and then set ad_option
% accordingly. On a 64-bit MATLAB, maxArrayDim = 2^48-1 according to the document of MATLAB R2019a.
% !!! Make sure that everything is compiled with the SAME ad_option !!!
% !!! Otherwise, Segmentation Fault may occur !!!
[Architecture, maxArrayDim] = computer;
if any(strfind(Architecture, '64')) && log2(maxArrayDim) > 31
    ad_option = '-largeArrayDims';
else
    ad_option = '-compatibleArrayDims'; % This works for both 32-bit and 64-bit MATLAB
end

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
    if debug_flags{idbg}
        mex_options = {ad_option, '-silent', '-g'};
    else
        mex_options = {ad_option, '-silent', '-O'};
    end
    for iprc = 1 : length(precisions)
        prepare_header(header_file, precisions{iprc}, debug_flags{idbg});
        work_dir = fullfile(common, pdstr(precisions{iprc}, debug_flags{idbg}));
        prepare_work_dir(work_dir);
        cd(work_dir);
        % One may write the loop below as
        %%mex(mex_options{:}, '-c', common_files{:});
        % But it does not work for some variants of MATLAB. This may be because the compilation above does
        % not respect the order of common_files{:}, which is critical due to the dependence among modules.
        for icf = 1 : length(common_files)
            mex(mex_options{:}, '-c', common_files{icf});
            % The module/object files are dumped to the current directory, namely `work_dir`.
        end
        common_obj_files = list_obj_files(work_dir);
        % Compile `gethuge`. We only provide a non-debugging version of `gethuge`.
        if ~debug_flags{idbg}
            gateway = fullfile(gateways, 'gethuge.F');
            mexname = get_mexname('gethuge', precisions{iprc});
            mex(mex_options{:}, common_obj_files{:}, gateway, '-output', mexname, '-outdir', mexdir);
        end
    end
end
fprintf('Done.\n');


% Compile the solvers.
for isol = 1 : length(solvers)
    solver = solvers{isol};
    fprintf('Compiling %s ... ', solver);
    gateway = fullfile(gateways, [solver(1:end-1), '_mex.F']);
    for ivar = 1 : length(variants)
        if strcmp(variants{ivar}, 'classical')
            soldir = fullfile(classical_src, solver(1:end-1));
        else
            soldir = fullfile(modern_src, solver(1:end-1));
        end
        for idbg = 1 : length(debug_flags)
            if strcmp(variants{ivar}, 'classical') && debug_flags{idbg}
                % The support for the classical variant is limited. No debugging version.
                continue
            end
            if debug_flags{idbg}
                mex_options = {ad_option, '-silent', '-g'};
            else
                mex_options = {ad_option, '-silent', '-O'};
            end
            for iprc = 1 : length(precisions)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~ismember(solver, ready_solvers) && ~strcmp(precisions{iprc}, 'double')
                    continue  % To be removed
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                work_dir = fullfile(soldir, pdstr(precisions{iprc}, debug_flags{idbg}));
                prepare_work_dir(work_dir);
                common_dir = fullfile(common, pdstr(precisions{iprc}, debug_flags{idbg}));
                copyfiles(list_mod_files(common_dir), work_dir);
                cd(work_dir);
                src_files = list_files(soldir, filelist);
                for isf = 1 : length(src_files)
                    mex(mex_options{:}, '-c', src_files{isf});
                    % The module/object files are dumped to the current directory, namely `work_dir`.
                end
                obj_files = [list_obj_files(common_dir), list_obj_files(work_dir)];
                mexname = get_mexname(solver, precisions{iprc}, debug_flags{idbg}, variants{ivar});
                mex(mex_options{:}, obj_files{:}, gateway, '-output', mexname, '-outdir', mexdir);
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
    rep_str(header_file, '#define __REAL_PRECISION__ 64', '#define __REAL_PRECISION__ 32');
    rep_str(header_file, '#define __REAL_PRECISION__ 128', '#define __REAL_PRECISION__ 32');
    rep_str(header_file, '#define __QP_AVAILABLE__ 1', '#define __QP_AVAILABLE__ 0');
case {'q', 'quadruple'}
    rep_str(header_file, '#define __REAL_PRECISION__ 32', '#define __REAL_PRECISION__ 128');
    rep_str(header_file, '#define __REAL_PRECISION__ 64', '#define __REAL_PRECISION__ 128');
    rep_str(header_file, '#define __QP_AVAILABLE__ 0', '#define __QP_AVAILABLE__ 1');
otherwise
    rep_str(header_file, '#define __REAL_PRECISION__ 32', '#define __REAL_PRECISION__ 64');
    rep_str(header_file, '#define __REAL_PRECISION__ 128', '#define __REAL_PRECISION__ 64');
    rep_str(header_file, '#define __QP_AVAILABLE__ 1', '#define __QP_AVAILABLE__ 0');
end

if debug_flag
    rep_str(header_file, '#define __DEBUGGING__ 0', '#define __DEBUGGING__ 1');
else
    rep_str(header_file, '#define __DEBUGGING__ 1', '#define __DEBUGGING__ 0');
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
