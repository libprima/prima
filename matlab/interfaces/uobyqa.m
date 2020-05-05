function [x, fx, exitflag, output] = uobyqa(varargin)
%UOBYQA is a solver for solving the following unconstrained continuous
%   optimization problem without using derivatives:
%
%   minimize    fun(x).
%
%   In the backend, UOBYQA calls late Professor M.J.D. Powell's Fotran code 
%   with the same name. The algorithm is described in [M. J. D. Powell,
%   UOBYQA: unconstrained optimization by quadratic approximation. Math.
%   Program., 92(B):555--582, 2002].
%
%   1. Basic syntax
%   
%   The command
%
%   x = uobyqa(fun, x0) 
%
%   solves the problem formulated above, where
%   *** fun is the name or function handle of the objective function
%   *** x0 is the starting point; x0 CANNOT be omitted or set to []
%
%   2. Outputs
%
%   The function can also be called with more outputs, e.g.,
%
%   [x, fx, exitflag, output] = uobyqa(INPUTS)
%
%   *** x is the approximate solution to the optimization problem
%   *** fx is fun(x)
%   *** exitflag is an integer indicating why UOBYQA returns; the
%       possible values are 
%       0: the lower bound for the trust region radius is reached
%       1: the target function value is achieved
%       2: a trust region step failed to reduce the quadratic model
%       3: the objective function has been evaluated maxfun times
%       4, 7, 8, 9: rounding errors become severe in the Fortran code 
%       -1: NaN occurs in x
%       -2: the objective function returns and NaN or nearly infinite
%       value (only in the classical mode)
%       -3: NaN occurs in the models
%       exitflag = 5, 10, 11, 12 are possible exitflags of the Fortran
%       code but cannot be returned by UOBYQA 
%   *** output is a structure with the following fields:
%       funcCount: number of function evaluations
%       fhist: history of function values
%       solver: backend solver that does the computation, i.e., 'uobyqa'
%       message: return message
%       warnings: a cell array that records all the  warnings raised
%       during the computation
%   
%   3. Options 
%
%   The same as FMINCON, UOBYQA accepts options passed by a structure.
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
%       of the approximate solution obtained by UOBYQA; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** classical: a boolean value indicating whether to call the classical 
%       Powell code or not; default: false
%   *** quiet: a boolean value indicating whether to keep quiet or not;
%       default: true (if false UOBYQA will print the return message of
%       the Fortran code)
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunval: a boolean value indicating whether to verify the returned 
%       function value or not; default: false
%       (if true, UOBYQA will check whether the returned value of fx
%       matches fun(x) or not, which costs a function evaluation) 
%
%   For example, the following code 
%   
%   options = struct();
%   options.maxfun = 50;
%   x = uobyqa(@cos, -1, options);
%
%   solves 
%       min cos(x) 
%   starting from x0=2 with at most 50 function evaluations.
%
%   4. Problem defined by a structure
%
%   The same as FMINCON, a problem can be passed to UOBYQA by a structure
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
%   x = uobyqa(problem);
%   
%   solves 
%       min cos(x) 
%   starting from x0=-1 with at most 50 function evaluations.
%
%   See also PDFO, NEWUOA, BOBYQA, LINCOA, COBYLA.
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
% Remarks: 
% !!! TREAT probinfo and options AS READONLY VARIABLES AFTER PREPDFO !!!
% !!! DO NOT MODIFY THE INFORMATION IN probinfo OR options AFTER PREPDFO !!! 
%
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% uobyqa starts

callstack = dbstack;
funname = callstack(1).name; % Name of the current function 
if length(callstack) >= 2
    invoker = callstack(2).name; % Name of the function who calls this function 
else
    invoker = '';
end
internal_invokers = {'pdfo'}; % Invokers from this package; may have others in the future

% OUTPUT records the information that is produced by the solver and
% intended to pass to postpdfo.
% OUTPUT should contain at least the following fiels:
% x, fx, exitflag, funcCount, fhist, constrviolation, chist, warnings;
% For lincoa, it should also contain constr_modified; for nonlinearly 
% constrained solvers, it should also contain nlcineq and nlceq. 
output = struct();
% N.B.: DO NOT record anything in PROBINFO or OPTIONS. If the solver is 
% called by pdfo, then postpdfo will do nothing; the real postprocessing 
% will be done when pdfo calls postpdfo using the OUTPUT returned by solver
% together with the PROBINFO and OPTIONS in pdfo; that said, in such a senario, 
% the PROBINFO and OPTIONS of this solver will NOT be passed to the real 
% postprocessing.

warning('off', 'backtrace'); % Do not display the stack trace of a warning
output.warnings = {}; % A cell that records all the warnings

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
% Even if invoker='pdfo', we still need to call prepdfo, which will assign 
% values to fun, x0, ..., options.
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
maxfun = options.maxfun;
rhobeg = options.rhobeg;
rhoend = options.rhoend;
ftarget = options.ftarget;

% Check whether the problem is too large for the Fortran code
% In the mex gateway, a workspace of size 
% nw = (n*(42+n*(23+n*(8+n)))+max(2*n*n+4,18*n))/4 + 1 (see below).
% will be allocated, which is the largest memory allocated by
% UOBYQA. If the value assigned to nw is so large that overflow
% occurs, then there will be a Segmentation Fault!!! 
% The largest possible value of nw depends on the type of nw in the
% mex file, which is the default INTEGER type (~2E9 for integer*4, 
% and ~9E18 for integer*8). This imposes an upper limit on the size the
% size of problem solvable by this code. If nw is INTEGER*4, the largest
% value of n is ~300. UOBYQA is not designed for so large problems.
% Indeed, when n > 10, NEWUOA/BOBYQA/LINCOA can solve unconstrained
% problems much more efficiently. 
% In the following code, gethuge returns the largest possible value of
% the given data type in the mex environment.

n = length(x0);

if (n <= 1)
    wid = sprintf('%s:UnivariateProblem', funname);
    wmessage = sprintf('%s: a univariate problem received; %s may fail. Try other solvers.', funname, funname);
    warning(wid, '%s', wmessage);
    output.warnings = [output.warnings, wmessage];
end


% The largest integer in the mex functions; the factor 0.99 provides a buffer
maxint = floor(0.99*min([gethuge('integer'), gethuge('mwSize'), gethuge('mwIndex')]));
nw = (n*(42+n*(23+n*(8+n)))+max(2*n*n+4,18*n))/4 + 1;
if nw >= maxint 
    % nw would suffer from overflow in the Fortran code; exit immediately 
    % Public/normal error
    if strcmp(invoker, 'pdfo')
        error(sprintf('%s:ProblemTooLarge', invoker), '%s: problem too large for %s. Try other solvers.', invoker, funname);
    else
        error(sprintf('%s:ProblemTooLarge', funname), '%s: problem too large for %s. Try other solvers.', funname, funname);
    end
end
if maxfun > maxint
    % maxfun would suffer from overflow in the Fortran code 
    maxfun = maxint;
    % Obviously, nw >= (n+1)(n+2)/2. If nw < maxint, 
    % then maxint > (n+1)(n+2)/2, and hence maxfun > (n+1)(n+2)/2 is guaranteed
    wid = sprintf('%s:MaxfunTooLarge', funname);
    wmessage = sprintf('%s: maxfun exceeds the upper limit of Fortran integers; it is set to %d.', funname, maxfun);
    warning(wid, '%s', wmessage);
    output.warnings = [output.warnings, wmessage];
end

try 
% Call the Fortran code
% The mexified Fortran Function is a private function generating only private errors; however, public errors can occur due to, e.g., evalobj; error handeling needed 
    if options.classical
        [x, fx, exitflag, nf, fhist] = fuobyqa_classical(fun, x0, rhobeg, rhoend, maxfun, ftarget);
    else
        [x, fx, exitflag, nf, fhist] = fuobyqa(fun, x0, rhobeg, rhoend, maxfun, ftarget);
    end
    % Record the results of the solver in OUTPUT
    output.x = x;
    output.fx = fx;
    output.exitflag = exitflag;
    output.funcCount = nf;
    output.fhist = fhist;
    output.constrviolation = 0; % Unconstrained problem
    output.chist = [];

% Postprocess the result 
% postpdfo are private functions that may generate public errors; error-handeling needed
    [x, fx, exitflag, output] = postpdfo(probinfo, output);
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

% uobyqa ends
return
