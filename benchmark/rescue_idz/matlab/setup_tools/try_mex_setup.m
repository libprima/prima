function success = try_mex_setup(language, verbose)
%TRY_MEX_SETUP tries running MEX setup for compiling language.
% At return,
% success = 1 means MEX is correctly set up.
% success = 0 means MEX setup fails.
% success = -1 means MEX setup runs successfully, but either we cannot try MEX on the example
% file because such a file is not found, or the MEX file of the example file works but the result
% is incorrect.

verbose = (nargin >=2 && verbose);

% Return if MEX is already well configured. This is important, because MEX is usable if it was set
% up before, and because MEX setup may fail even if it succeeded before due to change of environment.
success = mex_well_configured(language);  % verbose = false
if success == 1
    return
end

orig_warning_state = warning;
warning('off','all'); % We do not want to see warnings

% Try `mex('-setup', language)`
mex_setup = -1;
exception = [];
try
    %[~, mex_setup] = evalc('mex(''-setup'', language)'); % Use evalc so that no output will be displayed
    mex_setup = mex('-setup', language); % mex -setup may be interactive. So it is not good to mute it completely!!!
catch exception
    % Do nothing
end

% If MEX setup fails, it is probably because of failing to find a supported compiler. See
% https://www.mathworks.com/support/requirements/supported-compilers.html. On Linux, it is easy to
% fix and the user should be capable of handling it. So we focus on macOS and Windows.
% For C/C++, it suffices to install "Xcode with Clang" on macOS and "Microsoft Visual C++"
% (Microsoft Visual Studio with the "Desktop development with C++" workload) on Windows.
% For Fortran, we need Xcode on macOS or Microsoft Visual Studio (maybe also with "Desktop
% development with C++") on Windows, and additionally the Intel Fortran compiler with the
% environment variables ONEAPI_ROOT or IFORT_COMPILERYEAR set accordingly. In the following,
% for Fortran, we set the environment variables ONEAPI_ROOT and IFORT_COMPILER18/19 ...
% This should make MEX setup work for Fortran on MATLAB 2020a or above if a **default** installation
% of Intel OneAPI (available for free) has been done and Xcode or Microsoft VS is correctly installed.
% Update 20230401: MATLAB R2022b/R2023a seem not checking ONEAPI_ROOT, and IFORT_COMPILER23 is needed.
if strcmpi(language, 'fortran') && (ismac || ispc) && (~isempty(exception) || mex_setup ~= 0)
    if ismac
        oneapi_root = '/opt/intel/oneapi/';
        compiler_dir = [oneapi_root, 'compiler/latest/mac/'];
    elseif ispc  % Windows
        oneapi_root = 'C:\Program Files (x86)\Intel\oneAPI\';
        compiler_dir = [oneapi_root, 'compiler\latest\windows\'];
    end

    % Set PATH.
    compiler_bin = fullfile(compiler_dir, 'bin');
    compiler_bin64 = fullfile(compiler_bin, 'intel64');  % Why not worry about 32-bit case? Since R2016a, MATLAB has been 64-bit only.
    setenv('PATH', [getenv('PATH'), pathsep, compiler_bin, pathsep, compiler_bin64]);  % Not needed for Windows as of 2023.

    % Set IFORT_COMPILER18, IFORT_COMPILER19, ..., IFORT_COMPILERCURRENT, ONEAPI_ROOT.
    first_year = 18;
    current_year = year(datetime()) - 2000;
    nyear = current_year - first_year + 1;
    envvars = cell(1, nyear + 1);
    envvars_save = cell(1, nyear + 1);
    isenvs = false(1, nyear + 1);

    for ienvvar = 1 : length(envvars) - 1
        envvars{ienvvar} = ['IFORT_COMPILER', int2str(first_year + ienvvar - 1)];
    end
    envvars{end} = 'ONEAPI_ROOT';

    for ienvvar = 1 : length(envvars)
        envvar = envvars{ienvvar};
        % Test whether the environment variable exists (isenv is available since R2022b).
        isenvs(ienvvar) = ~exist('isenv', 'builtin') || isenv(envvar);
        % Save the value of the environment variable; the value is empty in case of nonexistence.
        envvars_save{ienvvar} = getenv(envvar);
        % Set the environment variable.
        if strcmp(envvar, 'ONEAPI_ROOT')
            setenv(envvar, oneapi_root);
        else
            setenv(envvar, compiler_dir);
        end
    end

    % Try setting up MEX again.
    mex_setup = -1;
    exception = [];
    try
        %[~, mex_setup] = evalc('mex(''-setup'', language)'); % Use evalc so that no output will be displayed
        % mex -setup may be interactive. So it is not good to mute it completely!!!
        if verbose
            mex_setup = mex('-v', '-setup', language);
        else
            mex_setup = mex('-setup', language);
        end
    catch exception
        % Do nothing
    end

    % If the setup fails again, give up after restoring the environment variables.
    if ~isempty(exception) || mex_setup ~= 0
        for ienvvar = 1 : length(envvars)
            envvar = envvars{ienvvar};
            setenv(envvar, envvars_save{ienvvar});
            if exist('unsetenv', 'builtin') && ~isenvs(ienvvar)  % unsetenv is available since R2022b.
                unsetenv(envvar);
            end
        end
    end
end

if ~isempty(exception) || mex_setup ~= 0
    fprintf('\nYour MATLAB failed to run mex(''-setup'', ''%s'').', language);
    fprintf('\nTo see the detailed error message, execute the following command:\n');
    fprintf('\n  mex(''-v'', ''-setup'', ''%s'')\n', language);
    success = 0;
else
    % Try `mex(example_file)`
    success = mex_well_configured(language, true);  % verbose = true
end

% Restore the behavior of displaying warnings
warning(orig_warning_state);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function success = mex_well_configured(language, verbose)
%MEX_WELL_CONFIGURED verifies whether MEX is well configured by testing it on the official example.
% At return,
% success = 1 means MEX compiles the example successfully and the resultant MEX file works well.
% success = 0 means MEX cannot compile the example or the resultant MEX file does not work.
% success = -1 means either we cannot try MEX on the example file because such a file is not found,
% or the MEX file of the example file works but the result is incorrect.

verbose = (nargin >=2 && verbose);

success = 1;

orig_warning_state = warning;
warning('off','all'); % We do not want to see warnings

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

% Locate example_file, which is an example provided by MATLAB for trying MEX.
% NOTE: MATLAB MAY CHANGE THE LOCATION OF THIS FILE IN THE FUTURE.
switch lower(language)
case 'fortran'
    example_file_name = 'timestwo.F';
case {'c', 'c++', 'cpp'}
    example_file_name = 'timestwo.c';
otherwise
    error(sprintf('%s:UnsupportedLang', funname), '%s: Language ''%s'' is not supported by %s.', funname, language, funname);
end
example_file = fullfile(matlabroot, 'extern', 'examples', 'refbook', example_file_name);

% Check whether example_file exists
if ~exist(example_file, 'file')
    if verbose
        fprintf('\n');
        wid = sprintf('%s:ExampleFileNotExist', funname);
        warning('on', wid);
        warning(wid, 'We cannot find\n%s,\nwhich is a MATLAB built-in example for trying MEX on %s. It will be ignored.\n', example_file, language);
    end
    success = -1;
    return
end

% Try `mex(example_file)`
%!------------------------------------------------------------------------------------------------!%
% In general, we should clear a MEX function before compiling it. Otherwise, it may lead to a
% failure of even crash. See https://github.com/equipez/test_matlab/tree/master/crash
% Without the next line, `mex(example_file)` fails on Windows if we run this script two times.
clear('timestwo');
%!------------------------------------------------------------------------------------------------!%
temp_mexdir = tempdir();  % The directory to output the MEX file of `timestwo`.
mex_status = -1;
exception = [];
try
    [~, mex_status] = evalc('mex(example_file, ''-outdir'', temp_mexdir)'); % Use evalc so that no output will be displayed
catch exception
    % Do nothing
end

if ~isempty(exception) || mex_status ~= 0
    delete(fullfile(temp_mexdir, 'timestwo.*'));  % Remove the trash before returning
    if verbose
        fprintf('\nThe MEX of your MATLAB failed to compile\n%s,\nwhich is a MATLAB built-in example for trying MEX on %s.\n', example_file, language);
        fprintf('\nTo see the detailed error message, execute the following command:\n');
        fprintf('\n  mex(''-v'', fullfile(matlabroot, ''extern'', ''examples'', ''refbook'', ''%s''));\n', example_file_name);
    end
    success = 0;
    return
end

% Check whether the mexified example_file works
addpath(temp_mexdir);  % Make `timestwo` available on path
exception = [];
try
    [~, timestwo_out] = evalc('timestwo(1)'); % Try whether timestwo works correctly
catch exception
    % Do nothing
end

rmpath(temp_mexdir);  % Clean up the path before returning.
delete(fullfile(temp_mexdir, 'timestwo.*'));  % Remove the trash before returning

if ~isempty(exception)
    if verbose
        fprintf('\nThe MEX of your MATLAB compiled\n%s,\nbut the resultant MEX file does not work.\n', example_file);
    end
    success = 0;
elseif abs(timestwo_out - 2) >= 20*eps
    if verbose
        fprintf('\n');
        wid = sprintf('%s:ExampleFileWorksIncorrectly', funname);
        warning('on', wid);
        warning(wid, 'The MEX of your MATLAB compiled\n%s,\nbut the resultant MEX file returns %.16f when calculating 2 times 1.', example_file, timestwo_out);
    end
    success = -1;
end

% Restore the behavior of displaying warnings
warning(orig_warning_state);

return
