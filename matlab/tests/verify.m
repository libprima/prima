function verify(varargin)
%This function tests the modernized version of Powell's solver against Powell's version, verifying
% whether they produce the same results on CUTEst problems.
%
% Usage:
% verify(solver, dimrange) , where `solver` is the name of the solver to test, while `dimrange`
% is the vector [mindim, maxdim], or "small", or "big", or "large";
% vector(solver, problem), where `problem` i the name of a problem to test;
% vector(solve, {problem, ir}), where `ir` is the index of the random run in `isequiv`.
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: July 2020
%
% Last Modified: Monday, October 04, 2021 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse the inputs.
[solver, options] = parse_input(varargin);

% Mexify the solvers.
mex_solvers(solver);

% Tell MATLAB where to find CUTEST.
locate_cutest();

% Conduct the verification.
solvers = {[solver, 'n'], solver};
if (isequiv(solvers, options))
    fprintf('\n\nThe test succeeds!\n\n');
else
    error(sfprintf('\n\nThe test FAILS!\n\n'));
end
