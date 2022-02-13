function setup(varargin)
%SETUP compiles the package and try adding the package into the search path.
%
%   Let solvername be a string indicating a solver name, and options be
%   a structure indicating compilation options. Then setup can be called
%   in the following ways:
%
%   setup(solvername, options)  % Compile a solver with options
%   setup(solvername)  % Compile a solver
%   setup(options)  % Compile all the solvers with options
%   setup path  % Add the paths needed to use the package
%   setup clean  % Remove all the compiled MEX files
%   setup uninstall  % Uninstall the package
%
%   REMARKS:
%
%   1. Since MEX is the standard way of calling Fortran code in MATLAB, you
%   need to have MEX properly configured for compile Fortran before using
%   the package. It is out of the scope of this package to help the users
%   to configure MEX.
%
%   2. At the end of this script, we will try saving the path of this package
%   to the search path. This can be done only if you have the permission to
%   write the following path-defining file:
%
%   fullfile(matlabroot, 'toolbox', 'local', 'pathdef.m')
%   NOTE: MATLAB MAY CHANGE THE LOCATION OF THIS FILE IN THE FUTURE
%
%   Otherwise, you CAN still use the package, except that you need to run
%   the startup.m script in the current directory each time you start a new
%   MATLAB session that needs the package. startup.m will not re-compile
%   the package but only add it into the search path.
%
%   ***********************************************************************
%   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University.
%
%   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: public (can be called directly by users)
%
% Remarks
%
% 1. Remarks on the directory mexdir.
% Functions and MEX files in the directory mexdir are
% automatically available to functions in the directory interfaces, and
% to scripts called by the functions that reside in interfaces. They are
% not available to other functions/scripts unless mexdir is
% added to the search path.
%
% 2. Remarks on the 'files_with_wildcard' function.
% MATLAB R2015b does not handle wildcard (*) properly. For example, if
% we would like to removed all the .mod files under a directory specified
% by dirname, then the following would workd for MATLAB later than R2016a:
% delete(fullfile(dirname, '*.mod'));
% However, MATLAB R2015b would complain that it cannot find '*.mod'.
% The 'files_with_wildcard' function provides a workaround.
%
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup starts

% Check the version of MATLAB.
if verLessThan('matlab', '8.3') % MATLAB R2014a = MATLAB 8.3
    fprintf('\nSorry, this package does not support MATLAB R2013b or earlier releases.\n\n');
    return
end

% All solvers to compile.
all_solvers = {'cobyla', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa'};

% The stamp to be included in the path_string. Needed only if `savepath` fails.
path_string_stamp = 'PDFO';

% The full path of several directories.
cpwd = fileparts(mfilename('fullpath')); % Current directory
fsrc = fullfile(cpwd, 'fsrc'); % Directory of the Fortran source code
fsrc_interform = fullfile(cpwd, 'fsrc', '.interform'); % Directory of the intersection-form Fortran source code
fsrc_common_interform = fullfile(fsrc_interform, 'common'); % Directory of the common files
fsrc_classical = fullfile(cpwd, 'fsrc', 'classical'); % Directory of the classical Fortran source code
fsrc_classical_interform = fullfile(cpwd, 'fsrc', 'classical', '.interform'); % Directory of the intersection-form Fortran source code
matd = fullfile(cpwd, 'matlab'); % Matlab directory
gateways = fullfile(matd, 'mex_gateways'); % Directory of the MEX gateway files
gateways_interform = fullfile(gateways, '.interform');  % Directory of the intersection-form MEX gateway files
interfaces = fullfile(matd, 'interfaces'); % Directory of the interfaces
mexdir = fullfile(interfaces, 'private'); % The private subdirectory of the interfaces
tests = fullfile(matd, 'tests'); % Directory containing some tests
tools = fullfile(matd, 'tools'); % Directory containing some tools, e.g., interform.m

% `tools` contains some functions needed below.
addpath(tools);

% Parse the input.
[solver_list, options, action, wrong_input] = parse_input(varargin);

% Exit if wrong input detected. Error messages have been printed during the parsing.
if wrong_input
    rmpath(tools);
    return
end

% Remove the compiled MEX files if requested.
if strcmp(action, 'clean')
    clean_mex;
    rmpath(tools);
    return
end

% Uninstall the package if requested.
if strcmp(action, 'uninstall')
    uninstall_pdfo(path_string_stamp);
    rmpath(tools);
    return
end

% Add the path and return if requested.
if strcmp(action, 'path')
    try
        add_save_path(interfaces, path_string_stamp);
    catch exception
        rmpath(tools);
        rethrow(exception);
    end
    rmpath(tools);
    fprintf('\nPath added.\n\n')
    return
end

% If we arrive here, then the user requests us to compile the solvers.

% Decide whether to compile with -O (optimize, the default) or with -g (debug).
% N.B.: -O and -g may lead to (slightly) different behaviors of the mexified code. This was observed
% on 2021-09-09 in a test of NEWUOA on the AKIVA problem of CUTEst. It was because the mexified code
% produced different results when it was supposed to evaluate COS(0.59843577329095299_DP) amid OTHER
% CALCULATIONS: with -O, the result was 0.82621783366991353; with -g, it became 0.82621783366991364.
% Bizarrely, if we write a short Fortran program to evaluate only COS(0.59843577329095299_DP),
% then the result is always 0.82621783366991364, regardless of -O or -g. No idea why.
if isempty(options)
    options = struct();
end
opt_option = '-O';  % Optimize the object code; this is the default
debug_flag = (isfield(options, 'debug') && options.debug);
if debug_flag
    opt_option = '-g';  % Debug mode; -g disables MEX's behavior of optimizing built object code
end

% Detect whether we are running a 32-bit MATLAB, where maxArrayDim = 2^31-1, and then set ad_option
% accordingly. On a 64-bit MATLAB, maxArrayDim = 2^48-1 according to the document of MATLAB R2019a.
% !!! Make sure that everything is compiled with the SAME ad_option !!!
% !!! Otherwise, Segmentation Fault may occur !!!
[Architecture, maxArrayDim] = computer;
if any(strfind(Architecture, '64')) && log2(maxArrayDim) > 31
    ad_option = '-largeArrayDims';
else
    ad_option = '-compatibleArrayDims'; % This will also work in a 64-bit MATLAB
end

% Set MEX options.
mex_options = [{opt_option}, {ad_option}, '-silent'];

% Check whether MEX is properly configured.
fprintf('\nVerifying the set-up of MEX ... \n\n');
language = 'FORTRAN'; % Language to compile
mex_well_conf = mex_well_configured(language);
if mex_well_conf == 0
    fprintf('\nVerification FAILED.\n')
    fprintf('\nThe MEX of your MATLAB is not properly configured for compiling Fortran.');
    fprintf('\nPlease configure MEX before using this package. Try ''help mex'' for more information.\n\n');
    return
elseif mex_well_conf == -1
    fprintf('\nmex(''-setup'', ''%s'') runs successfully but we cannot verify that MEX works properly.', language);
    fprintf('\nWe will try to continue.\n\n');
else
    fprintf('\nMEX is correctly set up.\n\n');
end

% Generate the intersection-form Fortran source code
% We need to do this because MEX accepts only the (obselescent) fixed-form Fortran code on Windows.
% Intersection-form Fortran code can be compiled both as free form and as fixed form.
fprintf('Refactoring the Fortran code ... ');
interform(fsrc);
interform(fsrc_classical);
interform(gateways);
fprintf('Done.\n\n');

% Clean up the directories fsrc and gateways before compilation.
% This is important especially if there was previously another compilation with a different ad_option.
% Without cleaning-up, the MEX files may be linked with wrong .mod or .o files, which can lead to
% serious errors including Segmentation Fault!
dir_list = {fsrc_common_interform, gateways_interform, mexdir};
for idir = 1 : length(dir_list)
    cellfun(@(filename) delete(filename), list_modo_files(dir_list{idir}));
end

% Name of the file that contains the list of Fortran files. There should be such a file in each
% Fortran source code directory, and the list should indicate the dependence among the files.
filelist = 'ffiles.txt';

% Compilation starts
fprintf('Compilation starts. It may take some time ...\n');
% Change directory to mexdir. All the intermediate files produced by the compilation (e.g., .mod)
% will be dumped to this directory. They will be removed when the compilation finishes.
cd(mexdir);

exception = [];
try
    % Compilation of the common files. They are shared by all solvers. We compile them only once.
    % gateways_interform/debug.F contains debugging subroutines tailored for MEX.
    copyfile(fullfile(gateways_interform, 'debug.F'), fsrc_common_interform);
    % ppf.h contains preprocessing directives. Set __DEBUGGING__ according to debug_flag.
    header_file = fullfile(fsrc_common_interform, 'ppf.h');
    header_file_bak = fullfile(fsrc_common_interform, 'ppf.h.bak');
    copyfile(header_file, header_file_bak);
    if debug_flag
        rep_str(header_file, '#define __DEBUGGING__ 0', '#define __DEBUGGING__ 1');
    else
        rep_str(header_file, '#define __DEBUGGING__ 1', '#define __DEBUGGING__ 0');
    end
    % Common Fortran source files.
    common_files = regexp(fileread(fullfile(fsrc_common_interform, filelist)), '\n', 'split');
    common_files = strtrim(common_files(~cellfun(@isempty, common_files)));
    common_files = fullfile(fsrc_common_interform, common_files);
    common_files = [common_files, fullfile(gateways_interform, 'fmxapi.F'), fullfile(gateways_interform, 'cbfun.F')];
    % The loop below may be written in one line as follows:
    %mex(mex_options{:}, '-c', common_files{:}, '-outdir', fsrc_common_interform);
    % But it does not work for some versions of MATLAB. This may be because the compilation above does
    % not respect the order of common_files{:}, which is critical due to the dependence among modules.
    for icf = 1 : length(common_files)
        mex(mex_options{:}, '-c', common_files{icf}, '-outdir', fsrc_common_interform);
    end
    common_obj_files = list_obj_files(fsrc_common_interform);

    % Compilation of function gethuge
    gateway = fullfile(gateways_interform, 'gethuge.F');
    mexname = 'gethuge';
    mex(mex_options{:}, common_obj_files{:}, gateway, '-output', mexname, '-outdir', mexdir);

    version_list = {'m', 'c'};  % m - modernized, c - classical
    for isol = 1 : length(solver_list)

        solver = solver_list{isol};
        fprintf('Compiling %s ... ', solver);

        gateway = fullfile(gateways_interform, [solver, '_mex.F']);

        for iver = 1 : length(version_list)
            switch version_list{iver}
            case 'm'
                srcdir = fullfile(fsrc_interform, solver);
                mexname = ['f', solver, 'n'];
            case 'c'
                srcdir = fullfile(fsrc_classical_interform, solver);
                mexname = ['f', solver, 'n_classical'];
            end

            % Clean up the source file directory
            cellfun(@(filename) delete(filename), list_modo_files(srcdir));
            % Compile
            src_files = regexp(fileread(fullfile(srcdir, filelist)), '\n', 'split');
            src_files = strtrim(src_files(~cellfun(@isempty, src_files)));
            src_files = fullfile(srcdir, src_files);
            for isf = 1 : length(src_files)
                mex(mex_options{:}, '-c', src_files{isf}, '-outdir', srcdir);
            end
            obj_files = [common_obj_files, list_obj_files(srcdir)];
            mex(mex_options{:}, obj_files{:}, gateway, '-output', mexname, '-outdir', mexdir);
            % Clean up the source file directory
            cellfun(@(filename) delete(filename), list_modo_files(srcdir));
        end
        fprintf('Done.\n');
    end

    % Clean up fsrc_common_interform.
    cellfun(@(filename) delete(filename), list_modo_files(fsrc_common_interform));
    % Clean up mexdir.
    cellfun(@(filename) delete(filename), list_modo_files(mexdir));

    % Remove the intersection-form Fortran files unless we are debugging.
    if ~debug_flag
        rmdir(fsrc_interform, 's');
        rmdir(fsrc_classical_interform, 's');
        rmdir(gateways_interform, 's');
    end

catch exception
    % Do nothing for the moment.
end

if exist(header_file_bak, 'file')  % Restore header_file
    movefile(header_file_bak, header_file);
end
cd(cpwd); % Change directory back to cpwd

if ~isempty(exception)  % Rethrow any exception caught above.
    rethrow(exception);
end

% Compilation ends successfully if we arrive here.
fprintf('Package compiled successfully!\n');

% Add `interfaces` to the MATLAB path, and then try saving the path.
path_saved = add_save_path(interfaces, path_string_stamp);

fprintf('\nThe package is ready to use.\n');
fprintf('\nYou may now try ''help pdfo'' for information on the usage of the package.\n');

if isempty(setdiff(all_solvers, solver_list))
    addpath(tests);
    fprintf('\nYou may also run ''testpdfo'' to test the package on a few examples.\n\n');
else
    fprintf('\n');
end

if ~path_saved  %  `add_save_path` failed to save the path.
    add_path_string = sprintf('addpath(''%s'');', interfaces);
    fprintf('***** To use the pacakge in other MATLAB sessions, do one of the following. *****\n\n');
    fprintf('- run ''savepath'' right now if you have the permission to do so.\n\n');
    fprintf('- OR add the following line to your startup script\n');
    fprintf('  (see https://www.mathworks.com/help/matlab/ref/startup.html for information):\n\n');
    fprintf('  %s\n\n', add_path_string);
    fprintf('- OR come to the current directory and run ''setup path'' when you need the package.\n\n');
end

rmpath(tools);
% setup ends
return
