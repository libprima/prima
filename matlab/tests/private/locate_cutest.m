function cpaths = locate_cutest()
%This function tells MATLAB where to find CUTEst. The following lines should be written according to
% the installation of CUTEst on this machine.

homedir = fileparts(fileparts(userpath));
cdir_dft = fullfile(homedir, 'local', 'cutesif', 'cutest');

cdir = getenv('CUTEST');
if isempty(cdir)
    cdir = cdir_dft;
    setenv('CUTEST', cdir);  % This is needed by `cutestdir`, which will be called by `macup`.
end

cmtools = fullfile(fileparts(cdir), 'mtools', 'msrc');
cpaths = {cmtools};

for ip = 1 : length(cpaths)
    addpath(cpaths{ip});
end
