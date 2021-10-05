function locate_cutest()
%This function tells MATLAB where to find CUTEst. The following lines should be written according to
% the installation of CUTEst on this machine.
setenv('CUTEST', '~/local/cutesif/cutest');
setenv('MASTSIF', '~/local/cutesif/sif');
addpath('~/local/cutesif/mtools/msrc');
