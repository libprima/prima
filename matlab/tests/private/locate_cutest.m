function cpaths = locate_cutest()
%This function tells MATLAB where to find CUTEst. The following lines should be written according to
% the installation of CUTEst on this machine.

cdir_dft = fullfile(homedir(), 'local', 'cutesif', 'cutest');

%cdir = getenv('CUTEST');
%if isempty(cdir)
%    cdir = cdir_dft;
%    setenv('CUTEST', cdir);  % This is needed by `cutestdir`, which will be called by `macup`.
%end

%cmtools = fullfile(fileparts(cdir), 'mtools', 'msrc');
cmtools = fullfile(fileparts(cdir_dft), 'mtools', 'msrc');
cpaths = {cmtools};

old_dir = cd(fileparts(cmtools))
clear('setup');
setup
cd(old_dir)

for ip = 1 : length(cpaths)
    addpath(cpaths{ip});
end
