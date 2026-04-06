function cmpaths = locate_matcutest(directory)
%This function finds where MatCUTEst (https://github.com/matcutest/matcutest) is installed, adds the
% paths needed for using MatCUTEst, and returns these paths in a cell array.
% We search at most 10 levels below the given directory, whose default value is the home directory.
% N.B.: As of 202508, MatCUTEst supports only Linux.

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
    % Use `find` to locate the directory containing the signature path. If multiple such directories
    % exist, print a warning, list all of them, and use the first one. We only search at most 13
    % levels below the given directory to avoid long searching time.
    [status, cmtools] = system(['find ', directory, ' -maxdepth 13 -wholename ', name_str, ' -type d -print']);

    if status ~= 0 || isempty(cmtools)
        error('locate_matcutest:MatCUTEstNotFound', 'MatCUTEst is not found under %s.', directory);
    end

    % Remove the leading and trailing white-space characters, including '\n'. Without doing this,
    % the following `strsplit` will create an extra empty string at the end of the cell array.
    cmtools = strtrim(cmtools);
    cmtools = strsplit(cmtools, '\n');
    if length(cmtools) > 1
        warning('locate_matcutest:MultipleMatCUTEstFound', 'Multiple MatCUTEst installations are found under %s.', directory);
    end
    fprintf('\nThe MatCUTEst installations found under %s:\n\n', directory);
    for i = 1 : length(cmtools)
        fprintf('\t%s\n', cmtools{i});
    end
    cmtools = cmtools{1};  % Use the first one if multiple are found.
    fprintf('\nThe MatCUTEst to be used:\n\n\t%s\n\n', cmtools);
end

cmpaths = {cmtools};  % There may be other paths to include in the future.

for ip = 1 : length(cmpaths)
    addpath(cmpaths{ip});
end
