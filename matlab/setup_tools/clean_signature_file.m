function clean_signature_file(directory)
%CLEAN_SIGNATURE_FILE removes the signature file in `directory`.

signature_file_name = '.signature';
signature_file = fullfile(directory, signature_file_name);

if exist(signature_file, 'file')
    delete(signature_file);
end

% CLEAN_SIGNATURE_FILE ends
return
