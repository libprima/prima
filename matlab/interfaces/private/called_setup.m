function cstp = called_setup()
%CALLED_SETUP checks whether setup.m has been called.

cpwd = fileparts(mfilename('fullpath'));  % The directory where this script resides.
signature_file_name = '.signature';
signature_file = fullfile(cpwd, signature_file_name);

cstp = exist(signature_file, 'file');

return
