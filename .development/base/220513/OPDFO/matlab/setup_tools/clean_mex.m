function clean_mex(directory, is_verbose)
%CLEAN_MEX removes all the compiled MEX files in `directory`.

verbose = false;
if nargin >= 2
    if islogical(is_verbose)
        verbose = is_verbose;
    else
        verbose = ischarstr(is_verbose) && (strcmp(is_verbose, 'verbose') || strcmp(is_verbose, 'v'));
    end
end

if verbose
    fprintf('\nRemoving the compiled MEX files (if any) in the following directory:\n\n');
    fprintf('%s\n\n', directory);
end

% List  the compiled MEX files
mex_files = files_with_wildcard(directory, ['*.', mexext]);

% Unload the compiled MEX files from memory.
% !!! Without this, MATLAB may crash in some cases! See https://github.com/equipez/test_matlab/tree/master/crash
for imf = 1 : length(mex_files)
    % Note that `mex_files` contains full path, but we only need the `mexname`, without the path or
    % extension. Thus we need to use `fileparts` to retrieve the `mexname` first. Otherwise, `clear`
    % will do nothing (it does not raise a warning when requested to clear something nonexistent).
    [~, mexename] = fileparts(mex_files{imf});
    clear(mexename);
end

% Remove the compiled MEX files
cellfun(@(filename) delete(filename), mex_files);

if verbose
    fprintf('Done.\n\n');
end

return
