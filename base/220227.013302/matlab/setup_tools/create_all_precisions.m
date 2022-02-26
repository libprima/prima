function create_all_precisions(options_or_directory)
%CREATE_ALL_PRECISIONS creates `all_precisions.m` under the directory containing this script
% according to `options_or_directory`. `all_precisions.m` should return a cell array containing
% the names of all the precisions ('double', 'single', 'quadruple') available for the Fortran
% solvers in this package. It is created in the following way.
%
% 0. We assume that 'double' is always available.
%
% 1. If `options_or_directory` is a string, then it will be interpreted as the full path to
% a directory. The return of `all_precisions.m` will reflect the precisions available under this
% directory for ALL the solvers (i.e., the return of `all_solvers()`). The availability is
% determined by `isavailable(directory, precision)`.
%
% 2. If `options_or_directory` is a structure (or empty), then it will be interpreted as compilation
% options. The return of `all_precisions.m` will reflect the precisions available after the compilation.

% Default values for the availability of 'single' and 'quadruple'. They are used only if
% `options_or_directory` is a structure (i.e., it is indeed the compilation options).
single_precision = true;  
quadruple_precision = false;  

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

% Parse the inputs.

% Move the existing version of `all_precisions.m` to `all_precisions.bak`.
mfiledir = fileparts(mfilename('fullpath'));  % The directory where this file resides.
allprec = 'all_precisions';
allprec_file = fullfile(mfiledir, [allprec, '.m']);
allprec_file_bak = fullfile(mfiledir, [allprec, '.bak']);
if exist(allprec_file_bak, 'file')
    % Remove `allprec_file_bak` if it exists. This is necessary. Otherwise, `allprec_file` may be
    % restored incorrectly in case of failure in the sequel.
    delete(allprec_file_bak);
end
if exist(allprec_file, 'file')
    movefile(allprec_file, allprec_file_bak);
end

try  % We use `try ... catch ...` in order to restore `allprec_file` in case of any failure.
    if nargin ~= 1
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid number of inputs.', funname);
    elseif ischarstr(options_or_directory)

        directory = options_or_directory;
        single_precision = isavailable(directory, 'single');
        quadruple_precision = isavailable(directory, 'quadruple');

    elseif isa(options_or_directory, 'struct')

        options = options_or_directory;
        if isfield(options, 'single') && islogicalscalar(options.single)
            single_precision = options.single;
        end
        if isfield(options, 'quadruple') && islogicalscalar(options.quadruple)
            quadruple_precision = options.quadruple;
        end
    elseif ~isempty(options_or_directory)
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: invalid input received', funname);
    end

    % Decide the precision list string.
    precision_list_string = '''double''';
    if single_precision
        precision_list_string = [precision_list_string, ', ''single'''];
    end
    if quadruple_precision
        precision_list_string = [precision_list_string, ', ''quadruple'''];
    end
    precision_list_string = sprintf('precision_list = {%s};', precision_list_string);

    % Create the file.
    fid = fopen(allprec_file, 'w');  % Open/create file for writing. Discard existing contents.
    if fid == -1
        error('Cannot create file %s.', allprec_file);
    end
    fprintf(fid, 'function precision_list = %s()\n', allprec);
    fprintf(fid, '%%%s ', upper(allprec));
    fprintf(fid, ' returns a cell array containing the names of all the precisions available for the\n');
    fprintf(fid, '%% Fortran solvers in this package.\n');
    fprintf(fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid, '%% This file is created automatically by \n%% %s.m at %s.\n', mfilename, datestr(datetime(), 'yymmdd.HH:MM:SS'));
    fprintf(fid, '%% NEVER EDIT IT, OR THE EARTH WILL EXPLODE.\n');
    fprintf(fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
    fprintf(fid, '%s', precision_list_string);
    fprintf(fid, '\n\nreturn');
    fclose(fid);

catch exception
    if exist(allprec_file_bak, 'file')
        movefile(allprec_file_bak, allprec_file);  % Restore the backup of `allprec_file` upon failure.
    end
    rethrow(exception);
end

% If we arrive here, `allprec_file` has been created successfully. Remove the backup.
if exist(allprec_file_bak, 'file')
    delete(allprec_file_bak);
end

% CREATE_ALL_PRECISIONS ends
return
