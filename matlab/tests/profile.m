function profile(varargin)
%This function profiles the modernized version of Powell's solver against Powell's version.
%
% Usage: profile(solver, dimrange) , where `solver` is the name of the solver to test, while `dimrange`
% is the vector [mindim, maxdim], or "small", or "big", or "large".
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: October 2021
%
% Last Modified: Monday, October 04, 2021 PM09:19:19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[solver, options] = parse_input(varargin);

% Mexify the solvers.
mex_solvers(solver);

% Tell MATLAB where to find CUTEST.
locate_cutest();

% Profile the solvers.
solvers = {[solver, 'n'], solver};
tic;
perfdata(solvers, options);
toc;
