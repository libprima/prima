function cmpaths = locate_matcutest(directory)
%This function finds where MatCUTEst (https://github.com/equipez/matcutest) is installed, add the
% paths needed for using MatCUTEst, and return these paths in a cell array.
% We search at most 3 levels below the given directory, whose default value is $HOME (as of 202301,
% MatCUTEst supports only Linux).

if nargin < 1
    directory = getenv('HOME');
end

[~, cmtools] = system(['find ', directory, ' -maxdepth 6 -wholename "*/matcutest/mtools/src" -type d -print -quit']);

if isempty(cmtools)
    error('locate_matcutest:MatCUTEstNotFound', 'MatCUTEst is not found under %s.', directory);
end

cmtools = strtrim(cmtools);  % Remove the leading and trailing white-space characters, including '\n'.

cmpaths = {cmtools};  % There may be other paths to include in the future.

for ip = 1 : length(cmpaths)
    addpath(cmpaths{ip});
end
