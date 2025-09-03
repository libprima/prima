function [path, mexname] = official_mex_example(language)
% This function returns the path to the official MEX example and the name of the MEX file without extension.
% NOTE: MATLAB MAY CHANGE THE LOCATION OF THIS FILE IN THE FUTURE.

callstack = dbstack;
funname = callstack(1).name; % Name of the current function

switch lower(language)
case 'fortran'
    example_file_name = 'timestwo.F';
case {'c', 'c++', 'cpp'}
    example_file_name = 'timestwo.c';
otherwise
    eid = sprintf('%s:UnsupportedLang', funname);
    error(eid, '%s: Language ''%s'' is not supported by %s.', funname, language, funname);
end
path = fullfile(matlabroot, 'extern', 'examples', 'refbook', example_file_name);

if ~exist(path, 'file')
    eid = sprintf('%s:ExampleFileNotExist', funname);
    error(eid, 'We cannot find\n%s,\nwhich is supposed to be a MATLAB built-in example for trying MEX on %s.', path, language);
end

[~, mexname] = fileparts(path); % File name without extension
