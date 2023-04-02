function setpath(path_string)

path(path_string);
orig_warning_state = warning;
warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path.

try
    sys_pathdef = fullfile(matlabroot(), 'toolbox', 'local', 'pathdef.m');
    savepath(sys_pathdef);
catch
    % Do nothing.
end

warning(orig_warning_state); % Restore the behavior of displaying warnings.
