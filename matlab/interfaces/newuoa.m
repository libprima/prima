function [x, fx, exitflag, output] = newuoa(varargin)
%NEWUOA is a solver for solving the following unconstrained continuous
%   optimization problem without using derivatives:
%
%   minimize    fun(x).
%
%   In the backend, NEWUOA calls late Professor M.J.D. Powell's Fotran code 
%   with the same name. The algorithm is described in [M. J. D. Powell,
%   The NEWUOA software for unconstrained optimization without derivatives, 
%   In Large-Scale Nonlinear Optimization, eds. G. Di Pillo and M. Roma, 
%   pages 255--297, Springer, New York, US, 2006].
%
%   1. Basic syntax
%     
%   The command
%
%   x = newuoa(fun, x0) 
%
%   solves the problem formulated above, where
%   *** fun is the name or function handle of the objective function
%   *** x0 is the starting point; x0 CANNOT be omitted or set to []
%
%   2. Outputs
%
%   The function can also be called with more outputs, e.g.,
%
%   [x, fx, exitflag, output] = newuoa(INPUTS)
%
%   *** x is the approximate solution to the optimization pronblem
%   *** fx is fun(x)
%   *** exitflag is an integer indicating why NEWUOA returns; the
%       possible values are 
%       0: the lower bound for the trust region radius is reached
%       1: the target function value is achieved
%       2: a trust region step failed to reduce the quadratic model
%       3: the objective function has been evaluated maxfun times
%       4, 7, 8, 9: rounding errors become severe in the Fortran code 
%       -1: NaN occurs in x
%       -2: the objective function returns an NaN or nearly infinite
%       value (only in the classical mode)
%       -3: NaN occurs in the models
%       exitflag = 5, 10, 11, 12 are possible exitflags of the Fortran
%       code but cannot be returned by NEWUOA 
%   *** output is a structure with the following fields:
%       funcCount: number of function evaluations
%       fhist: history of function values
%       solver: backend solver that does the computation, i.e., 'newuoa'
%       message: return message
%       warnings: a cell array that records all the  warnings raised
%       during the computation
%   
%   3. Options 
%
%   The same as FMINCON, NEWUOA accepts options passed by a structure.
%   Such a structure should be passed as an additional input appended to
%   the end of the input list in the basic syntax.
%
%   The options include
%   *** maxfun: maximal number of function evaluations; default: 500*length(x0)
%   *** ftarget: target function value; default: -Inf
%   *** rhobeg: initial trust-region radius; typically, rhobeg should 
%       be about one tenth of the greatest expected change to a variable; 
%       rhobeg should be positive; default: 1
%   *** rhoend: final trust region radius; rhoend reflects the precision
%       of the approximate solution obtained by NEWUOA; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** npt: number of interpolation points for constructing a model
%       default: 2*length(x0)+1
%   *** classical: a boolean value indicating whether to call the classical 
%       Powell code or not; default: false
%   *** quiet: a boolean value indicating whether to keep quiet or not;
%       default: true (if false NEWUOA will print the return message of
%       the Fortran code)
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunval: a boolean value indicating whether to verify the returned 
%       function value or not; default: false
%       (if true, NEWUOA will check whether the returned value of fx
%       matches fun(x) or not, which costs a function evaluation) 
%
%   For example, the following code 
%   
%   options = struct();
%   options.maxfun = 50;
%   x = newuoa(@cos, -1, options);
%
%   solves 
%       min cos(x) 
%   starting from x0=2 with at most 50 function evaluations.
%
%   4. Problem defined by a structure
%
%   The same as FMINCON, a problem can be passed to NEWUOA by a structure
%   PROBLEM containing the following fields: 
%   PROBLEM.objective, PROBLEM.x0, PROBLEM.options, where
%   PROBLEM.objective is the function name or function handle of the
%   objective function (corresponding to the input 'fun' mentioned above), 
%   and all the other fields correspond to the inputs introduced above with
%   the same names. 
%
%   For example, the following code 
%
%   problem = struct();
%   problem.objective = @cos;
%   problem.x0 = -1;
%   problem.options.maxfun = 50;
%   x = newuoa(problem);
%   
%   solves 
%       min cos(x) 
%   starting from x0=-1 with at most 50 function evaluations.
%
%   See also PDFO, UOBYQA, BOBYQA, LINCOA, COBYLA.
%
%   See https://www.pdfo.net for more information.
%
%   ***********************************************************************
%   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk) 
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University.
%
%   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: public (can  be called directly by users)
% 
% Remarks: None
%
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% newuoa starts
callstack = dbstack;
funname = callstack(1).name; % Name of the current function 
if length(callstack) >= 2
    invoker = callstack(2).name; % Name of the function who calls this function 
else
    invoker = '';
end
internal_invokers = {'pdfo'}; % Invokers from this package; may have others in the future

warning('off', 'backtrace'); % Do not display the stack trace of a warning
output.warnings = {}; % A cell that records all the warnings
% Why do we record the warning message in output.warnings
% instead of probinfo.warnings? Because, if newuoa is called by
% pdfo, then probinfo will not be passed to postpdfo, and hence
% the warning message will be lost. To the contrary, output will
% be passed to postpdfo anyway. 

maxarg = 3; % Maximal number of inputs
nvararg = length(varargin); % Number of inputs

% Interpret the input.
% Expected inputs: [fun, x0, options], yet some of them may be omitted.
if (nvararg < 1)
    if ismember(invoker, internal_invokers) % Private/unexpected error
        error(sprintf('%s:TooFewInputs', funname), '%s: UNEXPECTED ERROR: at least 1 input.', funname);
    else % Public/normal error
        error(sprintf('%s:TooFewInputs', funname), '%s: at least 1 input.', funname);
    end
elseif (nvararg == 1)
    args = varargin; % If there is only 1 input, then it is a structure specifying the problem
elseif (nvararg >= 2 && nvararg <= maxarg) 
    varargin = [varargin, cell(1, maxarg-nvararg)]; % 'augment' the inputs to maxarg by adding []
    % cell(m,n) returns an mxn array of [] 
    args = [varargin(1:2), cell(1, 7), varargin(end)]; % args{:} (should have 10 entries) will be the inputs for prepdfo
else
    if ismember(invoker, internal_invokers) % Private/unexpected error
        error(sprintf('%s:TooManyInputs', funname), '%s: UNEXPECTED ERROR: at most %d inputs.', funname, maxarg);
    else % Public/normal error
        error(sprintf('%s:TooManyInputs', funname), '%s: at most %d inputs.', funname, maxarg);
    end
end

% Preprocess the input 
try % prepdfo is a private function that may generate public errors; error-handeling needed
    [fun, x0, ~, ~, ~, ~, ~, ~, ~, options, probinfo] = prepdfo(args{:}); 
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

% Extract the options
npt = options.npt;
maxfun = options.maxfun;
rhobeg = options.rhobeg;
rhoend = options.rhoend;
ftarget = options.ftarget;

% Check whether the problem is too large for the Fortran code
% In the mex gateway, a workspace of size 
% nw = (npt+13)*(npt+n)+3*n*(n+3)/2 + 1 (see below).
% will be allocated, which is the largest memory allocated by
% NEWUOA. If the value assigned to nw is so large that overflow
% occurs, then there will be a Segmentation Fault!!! 
% The largest possible value of nw depends on the type of nw in the
% mex file, which is the default INTEGER type (~2E9 for integer*4, 
% and ~9E18 for integer*8). This imposes an upper limit on the size the
% size of problem solvable by this code. If nw is INTEGER*4, assuming
% that npt=2n+1, the largest value of n is ~16000. NEWUOA is not
% designed for so large problems.
% In the following code, gethuge returns the largest possible value of
% the given data type in the mex environment.

% The largest integer in the mex functions; the factor 0.99 provides a buffer
maxint = floor(0.99*min([gethuge('integer'), gethuge('mwSize'), gethuge('mwIndex')]));
n = length(x0);
minnw = (n+15)*(2*n+2)+3*n*(n+3)/2 + 1;
% minnw is the smallest possible nw, i.e., nw with the smallest npt, i.e., npt=n+2 
if minnw >= maxint  
    % nw would suffer from overflow in the Fortran code; exit immediately 
    % Public/normal error
    if strcmp(invoker, 'pdfo')
        error(sprintf('%s:ProblemTooLarge', invoker), '%s: problem too large for %s. Try other solvers.', invoker, funname);
    else
        error(sprintf('%s:ProblemTooLarge', funname), '%s: problem too large for %s. Try other solvers.', funname, funname);
    end
end
maxnpt = max(n+2, floor(0.5*(-(n+13)+sqrt((n-13)^2+4*(maxint-3*n*(n+3)/2-1)))));
% maxnpt is the largest possible value of npt given that nw <= maxint
if npt > maxnpt
    npt = maxnpt;
    wid = sprintf('%s:NptTooLarge', funname);
    wmessage = sprintf('%s: npt is so large that it is unable to allocate the workspace; it is set to %d.', funname, npt);
    warning(wid, '%s', wmessage);
    output.warnings = [output.warnings, wmessage];
end
if maxfun > maxint
    % maxfun would suffer from overflow in the Fortran code 
    maxfun = maxint;
    wid = sprintf('%s:MaxfunTooLarge', funname);
    wmessage = sprintf('%s: maxfun exceeds the upper limit of Fortran integers; it is set to %d.', funname, maxfun);
    warning(wid, '%s', wmessage);
    output.warnings = [output.warnings, wmessage];
end

try 
% Call the Fortran code
% The mexified Fortran Function is a private function generating only private errors; however, public errors can occur due to, e.g., evalobj; error handeling needed 
    if options.classical
        [x, fx, exitflag, nf, fhist] = fnewuoa_classical(fun, x0, rhobeg, rhoend, maxfun, npt, ftarget);
    else
        [x, fx, exitflag, nf, fhist] = fnewuoa(fun, x0, rhobeg, rhoend, maxfun, npt, ftarget);
    end

% Postprocess the result 
% postpdfo is a private function that may generate public errors; error-handeling needed
    [x, fx, exitflag, output] = postpdfo(x, fx, exitflag, output, nf, fhist, 0, [], options, probinfo);
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

% newuoa ends
return
