function create_signature_file(directory)
%CREATE_SIGNATURE_FILE creates a signature file in `directory`.

signature = datestr(datetime(), 'yymmdd.HH:MM:SS');
signature_file_name = '.signature';
signature_file = fullfile(directory, signature_file_name);

fid = fopen(signature_file, 'w');  % Open/create file for writing. Discard existing contents.
if fid == -1
    error('Cannot to open or create file %s', signature_file);
end
fprintf(fid, '%s', signature);
fclose(fid);

% CREATE_SIGNATURE_FILE ends
return
