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

% The stamp to be included in the path_string. Needed only if `savepath` fails.
package_name = 'PDFO';

% The full paths to several directories needed for the setup.
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
tools = fullfile(matd, 'setup_tools'); % Directory containing some tools, e.g., interform.m

% `tools` contains some functions needed in the sequel.
addpath(tools);
% The following files are needed after the setup. Copy them to `mexdir`.
copyfile(fullfile(tools, 'all_solvers.m'), mexdir);
copyfile(fullfile(tools, 'all_precisions.m'), mexdir);
copyfile(fullfile(tools, 'get_mexname.m'), mexdir);

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
    uninstall_pdfo(package_name);
    rmpath(tools);
    return
end

% Add the path and return if requested.
if strcmp(action, 'path')
    try
        add_save_path(interfaces, package_name);
    catch exception
        rmpath(tools);
        rethrow(exception);
    end
    rmpath(tools);
    fprintf('\nPath added.\n\n')
    return
end

% If we arrive here, then the user requests us to compile the solvers.

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

% Compilation starts
fprintf('Compilation starts. It may take some time ...\n');
% Change directory to mexdir. All the intermediate files produced by the compilation (e.g., .mod)
% will be dumped to this directory. They will be removed when the compilation finishes.
cd(mexdir);
exception = [];
try
    compile(solver_list, mexdir, fsrc_interform, fsrc_classical_interform, fsrc_common_interform, gateways_interform, options);
catch exception
    % Do nothing for the moment.
end

%% Remove the intersection-form Fortran files unless we are debugging.
%if ~debug_flag
%    rmdir(fsrc_interform, 's');
%    rmdir(fsrc_classical_interform, 's');
%    rmdir(gateways_interform, 's');
%end
cd(cpwd); % Change directory back to cpwd

if ~isempty(exception)
    rethrow(exception);  % Rethrow any exception caught during the compilation.
end

% Compilation ends successfully if we arrive here.
fprintf('Package compiled successfully!\n');

% Add `interfaces` to the MATLAB path, and then try saving the path.
path_saved = add_save_path(interfaces, package_name);

fprintf('\nThe package is ready to use.\n');
fprintf('\nYou may now try ''help pdfo'' for information on the usage of the package.\n');

if isempty(setdiff(all_solvers(), solver_list))
    addpath(tests);
    fprintf('\nYou may also run ''testpdfo'' to test the package on a few examples.\n\n');
else
    fprintf('\n');
end

if ~path_saved  %  `add_save_path` failed to save the path.
    add_path_string = sprintf('addpath(''%s'');', interfaces);
    fprintf('***** To use the package in other MATLAB sessions, do ONE of the following. *****\n\n');
    fprintf('- run ''savepath'' right now if you have the permission to do so.\n\n');
    fprintf('- OR add the following line to your startup script\n');
    fprintf('  (see https://www.mathworks.com/help/matlab/ref/startup.html for information):\n\n');
    fprintf('  %s\n\n', add_path_string);
    fprintf('- OR come to the current directory and run ''setup path'' when you need the package.\n\n');
end

rmpath(tools);

% setup ends
return
