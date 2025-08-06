function cfiles = copy_crash_dump_files(directory, verbose)
%This function copies the crash dump files to a specified directory.

if nargin < 2
    verbose = false;  % Default to not displaying the content of the crash dump files
end

pc = parcluster();
pc_dir = pc.JobStorageLocation;

% Below are possible locations of the crash dump files. The last three are for Windows.
files = dir(fullfile(pc_dir, '*', 'matlab_crash_dump*'));
files = [files; dir(fullfile(homedir(), 'matlab_crash_dump*'))];
files = [files; dir(fullfile(tempdir(), 'matlab_crash_dump*'))];
files = [files; dir(fullfile(getenv('temp'), 'matlab_crash_dump*'))];
files = [files; dir(fullfile('C:\Users\*\AppData\Local\Temp', 'matlab_crash_dump*'))];

cfiles = cell(length(files), 1);
for ifiles = 1 : length(files)
    cfiles{ifiles} = fullfile(files(ifiles).folder, files(ifiles).name);
    copyfile(cfiles{ifiles}, directory);
    if verbose
        type(cfiles{ifiles});  % Display the content of the crash dump file
    end
end

return
