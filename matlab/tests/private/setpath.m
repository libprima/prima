function setpath(path_string)

path(path_string);

try
    sys_pathdef = fullfile(matlabroot(), 'toolbox', 'local', 'pathdef.m');
    savepath(sys_pathdef);
catch
    % Do nothing.
end
