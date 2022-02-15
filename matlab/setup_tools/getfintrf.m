function succ = getfintrf(directory)
%GETFINTRF gets the `fintrf.h` file of the current MATLAB and save it to `directory`.

% "directory" can be given by a full path or a path relative to the
% current directory. The following lines get its full path.
if nargin < 1
    directory = cd();  % When "directory" is not given, we default it to the current directory
end
origdir = cd();
cd(directory);
directory = cd();  % Full path of the given directory, which is the current directory now.
cd(origdir);

% Copy `fintrf.h` to `directory`.
fintrf = fullfile(directory, 'fintrf.h');
copyfile(fullfile(matlabroot, 'extern', 'include', 'fintrf.h'), fintrf);
fileattrib(fintrf, '+w')

% Add time stamp to fintrf.
time = datestr(datetime(), 'HH.MM.SS, yyyy-mm-dd');
matv = [version, ', ', computer];
S = fileread(fintrf);
S = ['/* MATLAB version: ', matv, ' */', newline, '/* Retrieved at ', time, ' */', newline, S];
fid = fopen(fintrf, 'w');  % Open/create file for writing. Discard existing contents.
if fid == -1
    error('Cannot open file %s', fintrf);
end
fwrite(fid, S, 'char');
fclose(fid);
succ = true;

return
