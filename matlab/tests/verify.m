function verify(varargin)
%This function tests the modernized version of Powell's solver against Powell's version, verifying
% whether they produce the same results on CUTEst problems.
%
% Usage:
%
%   verify(solver, dimrange, nocompile_flag, sequential_flag, reverse_flag, problem_type, options)
%   verify(solver, problem, nocompile_flag, reverse_flag, problem_type, options)
%   verify(solve, problem, ir, nocompile_flag, reverse_flag, problem_type, options)
%
% where
% - `solver` is the name of the solver to test
% - `dimrange` is the vector [mindim, maxdim], or "small", or "big", or "large"
% - `problem` is the name of the problem to test
% - `ir` is the index of the random run in `isequiv`.
% - `nocompile_flag` is either 'nocompile' or 'ncp', indicating not to compile the solves
% - `sequential_flag` (optional) is either 'sequential' or 'seq', which means to test the problems sequentially
% - `reverse_flag` (optional) is either 'reverse' or 'rev', which means to test the solvers in the reverse order
% - `problem_type` can be any of {'u', 'b', 'l', 'n', 'ub', 'ubl', 'ubln', 'bl', 'bln', 'ln'},
%   indicating the problem type to test
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: July 2020
%
% Last Modified: Monday, February 12, 2022 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldpath = path();  % Record the current path.
restoredefaultpath;  % Restore the "right out of the box" path of MATLAB

olddir = pwd();  % Record the current directory.

% Parse the inputs.
[solver, options] = parse_input(varargin);


% Prepare the test directory, i.e., `test_dir`.
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
test_dir = prepare_test_dir(solver, funname, options);

exception = [];

try

    % Make the solvers available. Note that the solvers are under `test_dir`.
    get_solvers(solver, test_dir, options.compile);

    % Tell MATLAB where to find CUTEST.
    locate_cutest();

    % Go to the test directory. This is not really necessary. It will not affect the test, but any
    % output (e.g., NEWUOA_output.txt, fort.6) will be dumped to `test_dir`.
    cd(test_dir);

    % Record `olddir` in `options` so that we can come back to `olddir` during `isequiv` if
    % necessary (for example, when a single test fails).
    options.olddir = olddir;

    % Define the solvers to test.
    if isfield(options, 'reverse') && options.reverse
        solvers = {solver, [solver, 'n']};  % Reverse order: first run 'SOLVER', and then run 'SOLVERn'
    else
        solvers = {[solver, 'n'], solver};  % Default order: first run 'SOLVERn', and then run 'SOLVER'
    end

    % Show current path information.
    showpath(solvers);

    % Conduct the verification.
    tic; isequiv(solvers, options); toc;  % `isequiv` raises an error in case the solver behave differently.

    % Show current path information again at the end of test.
    showpath(solvers);

catch exception

    % Do nothing for the moment.

end

setpath(oldpath);  % Restore the path to oldpath.
cd(olddir);  % Go back to olddir.
fprintf('\nCurrently in %s\n\n', pwd());

if ~isempty(exception)  % Rethrow any exception caught above.
    rethrow(exception);
end
