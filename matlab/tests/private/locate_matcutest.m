function cmpaths = locate_matcutest(directory)
%This function finds where MatCUTEst (https://github.com/equipez/matcutest) is installed, adds the
% paths needed for using MatCUTEst, and returns these paths in a cell array.
% We search at most 3 levels below the given directory, whose default value is the home directory.
% N.B.: As of 202301, MatCUTEst supports only Linux.

% We use the following path as the signature to identify MatCUTEst.
signature_path = fullfile('matcutest', 'mtools', 'src');

% cmtools is the path to the directory containing the signature path.
cmtools = '';

if nargin < 1
    path_strs = strsplit(path(), pathsep);
    ind = find(endsWith(path_strs, signature_path), 1, 'first');
    if ~isempty(ind)
        cmtools = path_strs{ind};
    else
        directory = getenv('HOME');
    end
end

if isempty(cmtools)
    % In the following line, the "*/" before signature_path cannot be removed.
    name_str = ['"*/', signature_path, '"'];
    [~, cmtools] = system(['find ', directory, ' -maxdepth 6 -wholename ', name_str, ' -type d -print -quit']);

    if isempty(cmtools)
        error('locate_matcutest:MatCUTEstNotFound', 'MatCUTEst is not found under %s.', directory);
    end
end

cmtools = strtrim(cmtools);  % Remove the leading and trailing white-space characters, including '\n'.

cmpaths = {cmtools};  % There may be other paths to include in the future.

for ip = 1 : length(cmpaths)
    addpath(cmpaths{ip});
end
