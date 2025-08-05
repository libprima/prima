function path_saved = add_save_path(path_string, path_string_stamp)
%ADD_SAVE_PATH adds the path indicated by PATH_STRING to the MATLAB path and then tries saving path.
% PATH_STRING_STAMP is a stamp used when writing PATH_STRING to the user's startup.m file, which is
% needed only if `savepath` fails.

if nargin < 2
    path_string_stamp = sprintf('Added by %s', mfilename);
end

if ~exist(path_string, 'dir')
    error('prima_norma:PathNotExist', 'The string %s does not correspond to an existing directory.', path_string);
end

addpath(path_string);

% Try saving the path in the system path-defining file at sys_pathdef. If the user does not have
% writing permission for this file, then the path will not saved.
% N.B. Do not save the path to the pathdef.m file under userpath. This file is not loaded by default
% at startup. See
% https://www.mathworks.com/matlabcentral/answers/269482-is-userdata-pathdef-m-for-local-path-additions-supported-on-linux
orig_warning_state = warning;
warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path
sys_pathdef = fullfile(matlabroot(), 'toolbox', 'local', 'pathdef.m');
path_saved = (savepath(sys_pathdef) == 0);
warning(orig_warning_state); % Restore the behavior of displaying warnings

% If path not saved, try editing the startup.m of this user. Do this only if userpath is nonempty.
% Otherwise, we will only get a startup.m in the current directory, which will not be executed
% when MATLAB starts from other directories. On Linux, the default value of userpath is
% ~/Documents/MATLAB, but it will be '' if this directory does not exist. We refrain from creating
% this directory in that case.
if ~path_saved && numel(userpath) > 0
    user_startup = fullfile(userpath, 'startup.m');
    add_path_string = sprintf('addpath(''%s'');', path_string);
    full_add_path_string = sprintf('%s  %s %s', add_path_string, '%', path_string_stamp);

    % First, check whether full_add_path_string already exists in user_startup or not.
    if exist(user_startup, 'file')
        startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
        path_saved = any(strcmp(startup_text_cells, full_add_path_string));
    end

    if ~path_saved
        % We first check whether the last line of the user startup script is an empty line (or the
        % file is empty or even does not exist at all). If yes, we do not need to put a line break
        % before the path adding string.
        if exist(user_startup, 'file')
            startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
            last_line_empty = isempty(startup_text_cells) || (isempty(startup_text_cells{end}) && ...
                isempty(startup_text_cells{max(1, end-1)}));
        else
            last_line_empty = true;
        end

        file_id = fopen(user_startup, 'a');  % Open/create file for writing. Append data to the end.
        if file_id ~= -1 % If FOPEN cannot open the file, it returns -1; We keep silent if it fails.
            if ~last_line_empty  % The last line of user_startup is not empty
                fprintf(file_id, '\n');  % Add a new empty line
            end
            fprintf(file_id, '%s', full_add_path_string);
            fclose(file_id);
            % Check that full_add_path_string is indeed added to user_startup.
            if exist(user_startup, 'file')
                startup_text_cells = regexp(fileread(user_startup), '\n', 'split');
                path_saved = any(strcmp(startup_text_cells, full_add_path_string));
            end
        end
    end
end
