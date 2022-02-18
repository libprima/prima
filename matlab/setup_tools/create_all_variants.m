function create_all_variants(options_or_directory)
%CREATE_ALL_VARIANTS creates `all_variants.m` under the directory containing this script
% according to `options_or_directory`. `all_variants.m` should return a cell array containing
% the names of all the variants ('modern', 'classical') available for the Fortran solvers in
% this package. It is created in the following way.
%
% 0. We assume that 'modern' is always available.
%
% 1. If `options_or_directory` is a string, then it will be interpreted as the full path to
% a directory. The return of `all_variants.m` will reflect the variants available under this
% directory for ALL the solvers (i.e., the return of `all_solvers()`). The availability is
% determined by `isavailable(directory, variant)`.
%
% 2. If `options_or_directory` is a structure (or empty), then it will be interpreted as compilation
% options. The return of `all_variants.m` will reflect the variants available after the compilation.

classical_variant = true;  % Default value for the availability of `classical`.

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

% Parse the inputs.

% Move the existing version of `all_variants.m` to `all_variants.bak`. This MUST be done first,
% because `isavailable` will call `all_variants()` if it exists, which would lead to incorrect
% output of  `isavailable`.
cpwd = fileparts(mfilename('fullpath'));  % The directory where this file resides.
allvar = 'all_variants';
allvar_file = fullfile(cpwd, [allvar, '.m']);
allvar_file_bak = fullfile(cpwd, [allvar, '.bak']);
if exist(allvar_file_bak, 'file')
    % Remove `allvar_file_bak` if it exists. This is necessary. Otherwise, `allvar_file` may be
    % restored incorrectly in case of failure in the sequel.
    delete(allvar_file_bak);
end
if exist(allvar_file, 'file')
    movefile(allvar_file, allvar_file_bak);
end

try  % We use `try ... catch ...` in order to restore `allvar_file` in case of any failure.
    if nargin ~= 1
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid number of inputs.', funname);
    elseif ischarstr(options_or_directory)

        directory = options_or_directory;
        classical_variant = isavailable(directory, 'classical');

    elseif isa(options_or_directory, 'struct')

        options = options_or_directory;
        if isfield(options, 'classical') && islogicalscalar(options.classical)
            classical_variant = options.classical;
        end

    elseif ~isempty(options_or_directory)
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid input received', funname);
    end

    % Decide the variant list string.
    variant_list_string = '''modern''';
    if classical_variant
        variant_list_string = [variant_list_string, ', ''classical'''];
    end
    variant_list_string = sprintf('variant_list = {%s};', variant_list_string);

    % Create the file.
    fid = fopen(allvar_file, 'w');  % Open/create file for writing. Discard existing contents.
    if fid == -1
        error('Cannot create file %s.', allvar_file);
    end
    fprintf(fid, 'function variant_list = %s()\n', allvar);
    fprintf(fid, '%%%s ', upper(allvar));
    fprintf(fid, ' returns a cell array containing the names of all the variants available for the\n');
    fprintf(fid, '%% Fortran solvers in this package.\n');
    fprintf(fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid, '%% This file is created automatically by \n%% %s.m at %s.\n', mfilename, datestr(datetime(), 'yymmdd.HH:MM:SS'));
    fprintf(fid, '%% NEVER EDIT IT, OR THE EARTH WILL EXPLODE.\n');
    fprintf(fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
    fprintf(fid, '%s', variant_list_string);
    fprintf(fid, '\n\nreturn');
    fclose(fid);

catch exception
    if exist(allvar_file_bak, 'file')
        movefile(allvar_file_bak, allvar_file);  % Restore the backup of `allvar_file` upon failure.
    end
    rethrow(exception);
end

% If we arrive here, `allvar_file` has been created successfully. Remove the backup.
if exist(allvar_file_bak, 'file')
    delete(allvar_file_bak);
end

% CREATE_ALL_VARIANTS ends
return
