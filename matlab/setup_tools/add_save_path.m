function path_saved = add_save_path(path_string, path_string_stamp)
%ADD_SAVE_PATH adds the path indicated by PATH_STRING to the MATLAB path and then tries saving path.

if nargin < 2
    path_string_stamp = sprintf('Added by %s', mfilename);
end

addpath(path_string);

% Try saving path.
% SAVEPATH saves the current MATLABPATH in the path-defining file, which is located at:
% fullfile(matlabroot, 'toolbox', 'local', 'pathdef.m')
% It returns 0 if the file was saved successfully; 1 otherwise.
% If savepath fails (probably because we do not have the permission to write the above pathdef.m
% file), then we try saving the path to the user-specific pathdef.m file located at userpath.
% On linux, userpath = '$HOME/Documents/MATLAB'. However, if $HOME/Documents does not exist,
% then userpath = []. In this case, we will not save path to the user-specific pathdef.m file.
% Otherwise, we will only get a pathdef.m in the current directory, which will not be executed
% when MATLAB starts from other directories.
orig_warning_state = warning;
warning('off', 'MATLAB:SavePath:PathNotSaved'); % Maybe we do not have the permission to save path.
path_saved = (savepath == 0 || (numel(userpath) > 0 && savepath(fullfile(userpath, 'pathdef.m')) == 0));
warning(orig_warning_state); % Restore the behavior of displaying warnings

% If path not saved, try editing the startup.m of this user. Do this only if userpath is nonempty.
if ~path_saved && numel(userpath) > 0
    user_startup = fullfile(userpath, 'startup.m');
    add_path_string = sprintf('addpath(''%s'');', path_string);
    full_add_path_string = sprintf('%s\t%s %s', add_path_string, '%', path_string_stamp);

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
