function clean_mex(directory)
%CLEAN_MEX removes the compiled MEX files in `directory`.

fprintf('\nRemoving the compiled MEX files (if any) in the following directory:\n\n');
fprintf('%s\n\n', directory);

% Remove the compiled MEX files
mex_files = files_with_wildcard(directory, '*.mex*');
cellfun(@(filename) delete(filename), mex_files);

fprintf('Done.\n\n');
return
