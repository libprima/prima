function cpath = locate_cutest()
%This function tells MATLAB where to find CUTEst. The following lines should be written according to
% the installation of CUTEst on this machine.
setenv('CUTEST', '~/local/cutesif/cutest');
cpath = '~/local/cutesif/mtools/msrc';
addpath(cpath);
