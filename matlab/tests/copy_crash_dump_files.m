function cfiles = copy_crash_dump_files(directory)
%This function copies the crash dump files to a specified directory.

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
    copyfile(fullfile(files(ifiles).folder, files(ifiles).name), directory);
end

return
