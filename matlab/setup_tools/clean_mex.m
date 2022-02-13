function clean_mex
%CLEAN_MEX removes the compiled MEX files.

fprintf('\nRemoving the compiled MEX files (if any) ... ');
% The full path of several directories.
cpwd = fileparts(mfilename('fullpath')); % Current directory
matd = fileparts(cpwd); % Matlab directory
interfaces = fullfile(matd, 'interfaces'); % Directory of the interfaces
mexdir = fullfile(interfaces, 'private'); % The private subdirectory of the interfaces

% Remove the compiled MEX files
mex_files = files_with_wildcard(mexdir, '*.mex*');
cellfun(@(filename) delete(filename), mex_files);

fprintf('Done.\n\n');
return
