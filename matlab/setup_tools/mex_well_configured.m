function success = mex_well_configured(language)
%MEX_WELL_CONFIGURED verifies the set-up of MEX for compiling language

orig_warning_state = warning;
warning('off','all'); % We do not want to see warnings when verifying MEX

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

ulang = upper(language);

success = 1;
% At return,
% success = 1 means MEX is well configured,
% success = 0 means MEX is not well configured,
% success = -1 means "mex -setup" runs successfully, but either we cannot try
% it on the example file because such a file is not found, or the MEX file of
% the example file does not work as expected.

% Locate example_file, which is an example provided by MATLAB for trying MEX.
% NOTE: MATLAB MAY CHANGE THE LOCATION OF THIS FILE IN THE FUTURE.
switch ulang
case 'FORTRAN'
    example_file = fullfile(matlabroot, 'extern', 'examples', 'refbook', 'timestwo.F');
case {'C', 'C++', 'CPP'}
    example_file = fullfile(matlabroot, 'extern', 'examples', 'refbook', 'timestwo.c');
otherwise
    error(sprintf('%s:UnsupportedLang', funname), '%s: Language ''%s'' is not supported by %s.', funname, language, funname);
end

% Try `mex('-setup', ulang)`
mex_setup = -1;
exception = [];
try
    %[~, mex_setup] = evalc('mex(''-setup'', ulang)'); % Use evalc so that no output will be displayed
    mex_setup = mex('-setup', ulang); % mex -setup may be interactive. So it is not good to mute it completely!!!
catch exception
    % Do nothing
end
if ~isempty(exception) || mex_setup ~= 0
    fprintf('\nYour MATLAB failed to run mex(''-setup'', ''%s'').\n', language);
    success = 0;
    return
end

% Check whether example_file exists
if ~exist(example_file, 'file')
    fprintf('\n')
    wid = sprintf('%s:ExampleFileNotExist', funname);
    warning('on', wid);
    warning(wid, 'We cannot find\n%s,\nwhich is a MATLAB built-in example for trying MEX on %s. It will be ignored.\n', example_file, language);
    success = -1;
    return
end

% Try `mex(example_file)`
temp_mexdir = tempdir();
mex_status = -1;
exception = [];
try
%    [~, mex_status] = evalc('mex(example_file, ''-outdir'', temp_mexdir)'); % Use evalc so that no output will be displayed
    mex_status = mex(example_file, '-outdir', temp_mexdir)
catch exception
    % Do nothing
end

trash_files = files_with_wildcard(temp_mexdir, 'timestwo.*');

if ~isempty(exception) || mex_status ~= 0
    cellfun(@(filename) delete(filename), trash_files);  % Clean up the trash before returning
    fprintf('\nThe MEX of your MATLAB failed to compile\n%s,\nwhich is a MATLAB built-in example for trying MEX on %s.\n', example_file, language);
    success = 0;
    return
end

% Check whether the mexified example_file works
addpath(temp_mexdir);  % Make `timestwo` available on path
exception = [];
try
    %[~, timestwo_out] = evalc('timestwo(1)'); % Try whether timestwo works correctly
    timestwo_out = timestwo(1)
catch exception
    % Do nothing
end

rmpath(temp_mexdir);  % Clean up the path before returning.
cellfun(@(filename) delete(filename), trash_files);  % Clean up the trash before returning

if ~isempty(exception)
    fprintf('\nThe MEX of your MATLAB compiled\n%s,\nbut the resultant MEX file does not work.\n', example_file);
    success = 0;
elseif abs(timestwo_out - 2)/2 >= 10*eps
    fprintf('\n')
    wid = sprintf('%s:ExampleFileWorksIncorrectly', funname);
    warning('on', wid);
    warning(wid, 'The MEX of your MATLAB compiled\n%s,\nbut the resultant MEX file returns %.16f when calculating 2 times 1.', example_file, timestwo_out);
    success = -1;
end

warning(orig_warning_state); % Restore the behavior of displaying warnings
return
