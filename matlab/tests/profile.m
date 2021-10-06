function profile(varargin)
%This function profiles the modernized version of Powell's solver against Powell's version.
%
% Usage:
%
%   profile(solver, dimrange, nocompile_flag)
%   profile(solver, dimrange, reload_flag)
%
% where
% - `solver` is the name of the solver to test
% - `dimrange` is the vector [mindim, maxdim], or "small", or "big", or "large"
% - `nocompile_flag` is either 'nocompile' or 'ncp', indicating not to compile the solves
% - `reload_flag` is either 'reload' or 'load', indicating not load the data directly from the .mat
% file corresponding to `solver` and `dimrange`
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: October 2021
%
% Last Modified: Monday, October 04, 2021 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse the inputs.
[solver, options] = parse_input(varargin);

% Mexify the solvers.
if options.compile
    mex_solvers(solver);
end

% Tell MATLAB where to find CUTEST.
locate_cutest();

% Set up the directory to save the testing data.
testdir = fileparts(mfilename('fullpath')); % Directory where this .m file resides.
datadir = fullfile(testdir, 'testdata');
if ~exist(datadir, 'dir')
    mkdir(datadir);
end
options.datadir = datadir;

% Profile the solvers.
solvers = {[solver, 'n'], solver};
tic;
perfdata(solvers, options);
toc;
