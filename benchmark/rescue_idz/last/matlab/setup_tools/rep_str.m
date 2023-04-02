function rep_str(filename, old_str, new_str)
%REP_STR replaces all `old_str` in filename with `new_str`.

fid = fopen(filename, 'r');  % Open file for reading.
if fid == -1
    error('Cannot open file %s.', filename);
end

% Read the file into a cell of strings
data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
fclose(fid);
cstr = data{1};

% Replace `old_str` with `new_str`.
for i = 1 : length(cstr)
    cstr{i} = strrep(cstr{i}, old_str, new_str);
end

% Save the file again
fid = fopen(filename, 'w');  % Open/create file for writing. Discard existing contents, if any.
if fid == -1
    error('Cannot open file %s.', filename);
end
fprintf(fid, '%s\n', cstr{:});
fclose(fid);

return
