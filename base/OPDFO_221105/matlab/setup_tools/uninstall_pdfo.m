function uninstall_pdfo(path_string_stamp)
%UNINSTALL_PDFO uninstalls PDFO.

fprintf('\nUninstalling PDFO (if it is installed) ... ');

% The full path of several directories.
mfiledir = fileparts(mfilename('fullpath'));  % The directory where this .m file resides
matd = fileparts(mfiledir); % Matlab directory
interfaces = fullfile(matd, 'interfaces'); % Directory of the interfaces
mexdir = fullfile(interfaces, 'private'); % The private subdirectory of the interfaces
tests = fullfile(matd, 'tests'); % Directory containing some tests

% Remove the compiled MEX files.
clean_mex(mexdir);

% Try removing the paths possibly added by PDFO
orig_warning_state = warning;
warning('off', 'MATLAB:rmpath:DirNotFound'); % Maybe the paths were not added. We do not want to see this warning.
warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path.
rmpath(interfaces, tests);
savepath;
warning(orig_warning_state); % Restore the behavior of displaying warnings

% Removing the line possibly added to the user startup script
user_startup = fullfile(userpath,'startup.m');
if exist(user_startup, 'file')
    add_path_string = sprintf('addpath(''%s'');', interfaces);
    full_add_path_string = sprintf('%s\t%s %s', add_path_string, '%', path_string_stamp);
    try
        del_str_ln(user_startup, full_add_path_string);
    catch
        % Do nothing.
    end
end

callstack = dbstack('-completenames');
root_dir = fileparts(callstack(2).file);  % Root directory of the package
fprintf('Done.\nYou may now remove\n\n    %s\n\nif it contains nothing you want to keep.\n\n', root_dir);

return
