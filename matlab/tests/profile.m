function profile(varargin)
%This function profiles the modernized version of Powell's solver against Powell's version.
%
% Usage:
%
%   profile(solver, dimrange, nocompile_flag, sequential_flag, reverse_flag, problem_type, options)
%   profile(solver, dimrange, reload_flag, reverse_flag, problem_type, options)
%
% where
% - `solver` is the name of the solver to test
% - `dimrange` is the vector [mindim, maxdim], or "small", or "big", or "large"
% - `nocompile_flag` is either 'nocompile' or 'ncp', indicating not to compile the solves
% - `reload_flag` is either 'reload' or 'load', indicating to load the data directly from the .mat
%   file corresponding to `solver` and `dimrange`
% - `sequential_flag` (optional) is either 'sequential' or 'seq', which means to test the problems sequentially
% - `reverse_flag` (optional) is either 'reverse' or 'rev', which means to test the solvers in the reverse order
% - `problem_type` can be any of {'u', 'b', 'l', 'n', 'ub', 'ubl', 'ubln', 'bl', 'bln', 'ln'},
%   indicating the problem type to test
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: October 2021
%
% Last Modified: Monday, October 04, 2021 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

olddir = cd();  % Record the current path.
oldpath = path();  % Record the current dir.
restoredefaultpath; % Restore the "right out of the box" path of MATLAB

try

    % Parse the inputs.
    [solver, options] = parse_input(varargin);

    % Make the solvers available.
    if options.compile
        % Mexify the solvers.
        mex_solvers(solver);
    else
        % Tell MATLAB where to find the solvers.
        locate_solvers();
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

    % Show current path information.
    showpath();

    % Profile the solvers.
    if isfield(options, 'reverse') && options.reverse
        solvers = {solver, [solver, 'n']};
    else
        solvers = {[solver, 'n'], solver};
    end
    tic;
    perfdata(solvers, options);
    toc;

catch exception

    setpath(oldpath);  % Restore the path to oldpath.
    cd(olddir);  % Go back to olddir.
    rethrow(exception);

end

setpath(oldpath);  % Restore the path to oldpath.
cd(olddir);  % Go back to olddir.
