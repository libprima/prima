function verify(varargin)
%This function tests the modernized version of Powell's solver against Powell's version, verifying
% whether they produce the same results on CUTEst problems.
%
% Usage:
%
%   verify(solver, dimrange, nocompile_flag, problem_type, options)
%   verify(solver, problem, nocompile_flag, problem_type, options)
%   verify(solve, problem, ir, nocompile_flag, problem_type, options)
%
% where
% - `solver` is the name of the solver to test
% - `dimrange` is the vector [mindim, maxdim], or "small", or "big", or "large"
% - `problem` is the name of the problem to test
% - `ir` is the index of the random run in `isequiv`.
% - `nocompile_flag` is either 'nocompile' or 'ncp', indicating not to compile the solves
% - `problem_type` can be any of {'u', 'b', 'l', 'n', 'ub', 'ubl', 'ubln', 'bl', 'bln', 'ln'},
%   indicating the problem type to test
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: July 2020
%
% Last Modified: Monday, October 04, 2021 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

olddir = cd();  % Record the current path.
oldpath = path();  % Record the current dir.
restoredefaultpath; % Restore the "right out of the box" path of MATLAB

try

    % Parse the inputs.
    [solver, options] = parse_input(varargin);

    % Mexify the solvers.
    if options.compile
        mex_solvers(solver);
    end

    % Tell MATLAB where to find CUTEST.
    locate_cutest();

    % Conduct the verification.
    solvers = {[solver, 'n'], solver};
    isequiv(solvers, options);  % `isequiv` raises an error in case the solver behave differently.

catch exception

    setpath(oldpath);  % Restore the path to oldpath.
    cd(olddir);  % Go back to olddir.
    rethrow(exception);

end

setpath(oldpath);  % Restore the path to oldpath.
cd(olddir);  % Go back to olddir.
