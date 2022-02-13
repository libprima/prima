function del_str_ln(filename, string)
%DEL_STR_LN deletes from filename all the lines that are identical to string
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file %s.', filename);
end

% Read the file into a cell of strings
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);
cstr = data{1};
% Remove the rows containing string
cstr(strcmp(cstr, string)) = [];

% Save the file again
fid = fopen(filename, 'w');
if fid == -1
    error('Cannot open file %s.', filename);
end
fprintf(fid, '%s\n', cstr{:});
fclose(fid);
return
