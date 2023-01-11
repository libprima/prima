function full_files = files_with_wildcard(dir_name, wildcard_string)
%FULL_FILES returns a cell array of files that match the wildcard_string under dir_name.
files = dir(fullfile(dir_name, wildcard_string));
full_files = cellfun(@(s)fullfile(dir_name, s), {files.name}, 'uniformoutput', false);
return
