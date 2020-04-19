function [x, fx, exitflag, output] = bobyqa(varargin)
%BOBYQA is a solver for solving the following bound-constrained continuous
%   optimization problem without using derivatives:
%
%   minimize    fun(x)
%       s.t.    lb <= x <= ub.
%
%   In the backend, BOBYQA calls late Professor M.J.D. Powell's Fotran code 
%   with the same name. The algorithm is described in [M. J. D. Powell,
%   The BOBYQA algorithm for bound constrained optimization without derivatives, 
%   Technical Report DAMTP 2009/NA06, Department of Applied Mathematics and 
%   Theoretical Physics, Cambridge University, Cambridge, UK, 2009].
%
%   1. Basic syntax
%
%   The command 
%
%   x = bobyqa(fun, x0, lb, ub) 
%
%   solves the problem formulated above, where
%   *** fun is the name or function handle of the objective function
%   *** x0 is the starting point; x0 CANNOT be []
%   *** lb and ub, which are vectors of the same length as x, are the 
%       lower and upper bound in the bound constraint lb <= x <= ub; 
%       set lb = [] if no lower bound, and ub = [] if no upper bound 
%
%   The function can also be called with more outputs, e.g.,
%
%   [x, fx, exitflag, output] = bobyqa(INPUTS)
%
%   See "3. Outputs" below for explanations on these outputs. 
%
%   2. Flexible syntax
%
%   x = bobyqa(fun, x0) solves 
%       minimize fun(x) 
%   x = bobyqa(fun, x0, lb) solves
%       minimize fun(x) s.t. lb <= x
%
%   3. Outputs
%
%   *** x is the approximate solution to the optimization pronblem
%   *** fx is fun(x)
%   *** exitflag is an integer indicating why BOBYQA returns; the
%       possible values are 
%       0: the lower bound for the trust region radius is reached
%       1: the target function value is achieved
%       2: a trust region step failed to reduce the quadratic model
%       3: the objective function has been evaluated maxfun times
%       4, 7, 8, 9: rounding errors become severe in the Fortran code 
%       13: all variables are fixed by the constraints
%       -1: NaN occurs in x
%       -2: the objective function returns an NaN or nearly infinite
%       value (only in the classical mode)
%       -3: NaN occurs in the models
%       -4: constraints are infeasible 
%       exitflag = 5, 10, 11, 12 are possible exitflags of the Fortran
%       code but cannot be returned by BOBYQA 
%   *** output is a structure with the following fields:
%       funcCount: number of function evaluations
%       constrviolation: constrviolation of x (if problem is
%       constrained; should be 0 since BOBYQA is a feasible method)
%       fhist: history of function values
%       chist: history of constraint violations (should be all 0)
%       solver: backend solver that does the computation, i.e., 'bobyqa'
%       message: return message
%       warnings: a cell array that records all the  warnings raised
%       during the computation
%   
%   4. Options 
%
%   The same as FMINCON, BOBYQA accepts options passed by a structure.
%   Such a structure should be passed as an additional input appended to
%   the end of the input list in the basic syntax or the flexible syntax. 
%
%   The options include
%   *** maxfun: maximal number of function evaluations; default: 500*length(x0)
%   *** ftarget: target function value; default: -Inf
%   *** rhobeg: initial trust-region radius; typically, rhobeg should 
%       be about one tenth of the greatest expected change to a variable; 
%       rhobeg should be positive; default: 1
%   *** rhoend: final trust region radius; rhoend reflects the precision
%       of the approximate solution obtained by BOBYQA; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** npt: number of interpolation points for constructing a model
%       default: 2*length(x0)+1
%   *** classical: a boolean value indicating whether to call the classical 
%       Powell code or not; default: false
%   *** scale: a boolean value that indicating whether to scale the problem
%       according to bounds or not; default: false
%   *** quiet: a boolean value indicating whether to keep quiet or not;
%       default: true (if false BOBYQA will print the return message of
%       the Fortran code)
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunval: a boolean value indicating whether to verify the returned 
%       function value or not; default: false
%       (if true, BOBYQA will check whether the returned value of fx
%       matches fun(x), which costs a function evaluation) 
%
%   For example, the following code 
%   
%   options = struct();
%   options.maxfun = 50;
%   x = bobyqa(@cos, -1, 2, 3, options);
%
%   solves 
%       min cos(x) s.t. 2 <= x <= 3 
%   starting from x0=2 with at most 50 function evaluations.
%
%   5. Problem defined by a structure
%
%   The same as FMINCON, a problem can be passed to BOBYQA by a structure
%   PROBLEM containing the following fields: 
%   PROBLEM.objective, PROBLEM.x0, PROBLEM.lb, PROBLEM.ub, PROBLEM.options, 
%   where PROBLEM.objective is the function name or function handle of
%   the objective function (corresponding to the input 'fun' mentioned above), 
%   and all the other fields correspond to the inputs introduced above with
%   the same names. 
%
%   For example, the following code 
%
%   problem = struct();
%   problem.objective = @cos;
%   problem.x0 = -1;
%   problem.lb = 2;
%   problem.ub = 3;
%   problem.options.maxfun = 50;
%   x = bobyqa(problem);
%   
%   solves 
%       min cos(x) s.t. 2 <= x <= 3 
%   starting from x0=-1 with at most 50 function evaluations.
%
%   See also PDFO, UOBYQA, NEWUOA, LINCOA, COBYLA.
%
%   See https://www.pdfo.co for more information.
%
%   ***********************************************************************
%   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk) 
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
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

% bobyqa starts

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
% instead of probinfo.warnings? Because, if bobyqa is called by
% pdfo, then probinfo will not be passed to postpdfo, and hence
% the warning message will be lost. To the contrary, output will
% be passed to postpdfo anyway. 

maxarg = 5; % Maximal number of inputs
nvararg = length(varargin); % Number of inputs

% Interpret the input.
% Expected inputs: [fun, x0, lb, ub, options], yet some of them may be omitted.
if (nvararg < 1)
    if ismember(invoker, internal_invokers) % Private/unexpected error
        error(sprintf('%s:TooFewInputs', funname), '%s: UNEXPECTED ERROR: at least 1 input.', funname);
    else % Public/normal error
        error(sprintf('%s:TooFewInputs', funname), '%s: at least 1 input.', funname);
    end
elseif (nvararg == 1)
    args = varargin; % If there is only 1 input, then it is a structure specifying the problem
elseif (nvararg >= 2 && nvararg <= maxarg) 
    % If 2<=nvararg<=5 and the last input is a structure (or []), then it is the 'options'
    if isa(varargin{end}, 'struct') 
        varargin = [varargin(1:end-1), cell(1, maxarg-nvararg), varargin(end)]; % 'augment' the inputs to maxarg by adding []
        % cell(m,n) returns an mxn array of [] 
    else
        varargin = [varargin, cell(1, maxarg-nvararg)]; % 'augment' the inputs to maxarg by adding []
    end
    args = [varargin(1:2), cell(1, 4), varargin(3:4), {[]}, varargin(end)]; % args{:} (should have 10 entries) will be the inputs for prepdfo
else
    if ismember(invoker, internal_invokers) % Private/unexpected error
        error(sprintf('%s:TooManyInputs', funname), '%s: UNEXPECTED ERROR: at most %d inputs.', funname, maxarg);
    else % Public/normal error
        error(sprintf('%s:TooManyInputs', funname), '%s: at most %d inputs.', funname, maxarg);
    end
end

% Preprocess the input 
try % prepdfo is a private function that may generate public errors; error-handeling needed
    [fun, x0, ~, ~, ~, ~, lb, ub, ~, options, probinfo] = prepdfo(args{:}); 
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

if probinfo.infeasible % The problem turned out infeasible during prepdfo
    exitflag = -4;
    nf = 0;
    x = NaN(size(x0));
    fx = NaN;
    fhist = [];
    constrviolation = NaN;
    chist = [];
elseif probinfo.nofreex % x was fixed by the bound constraints during prepdfo
    exitflag = 13;
    nf = 1;
    x = probinfo.fixedx_value;
    fx = fun(x);
    fhist = fx;
    constrviolation = probinfo.constrv_fixedx;
    chist = constrviolation;
else % The problem turns out 'normal' during prepdfo
    % Extract the options
    npt = options.npt;
    maxfun = options.maxfun;
    rhobeg = options.rhobeg;
    rhoend = options.rhoend;
    ftarget = options.ftarget;
    
    % Check whether the problem is too large for the Fortran code
    % In the mex gateway, a workspace of size 
    % nw = (npt+5)*(npt+n)+3*n*(n+5)/2 + 1
    % will be allocated, which is the largest memory allocated by
    % BOBYQA. If the value assigned to nw is so large that overflow
    % occurs, then there will be a Segmentation Fault!!! 
    % The largest possible value of nw depends on the type of nw in the
    % mex file, which is the default INTEGER type (~2E9 for integer*4, 
    % and ~9E18 for integer*8). This imposes an upper limit on the size
    % of problem solvable by this code. If nw is integer*4, assuming
    % that npt=2n+1, the largest value of n is ~16000. BOBYQA is not
    % designed for so large problems. 
    % In the following code, gethuge returns the largest possible value
    % of the given data type in the mex environment.
    
    % The largest integer in the mex functions; the factor 0.99 provides a buffer
    maxint = floor(0.99*min([gethuge('integer'), gethuge('mwSize'), gethuge('mwIndex')]));
    n = length(x0);
    minnw = (n+7)*(2*n+2)+3*n*(n+5)/2+1;
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
    maxnpt = max(n+2, floor(0.5*(-(n+5)+sqrt((n-5)^2+4*(maxint-3*n*(n+5)/2-1)))));
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

    % Call the Fortran code
    try % The mexified Fortran function is a private function generating only private errors; however, public errors can occur due to, e.g., evalobj; error handeling needed 
        if options.classical 
            [x, fx, exitflag, nf, fhist, constrviolation, chist] = fbobyqa_classical(fun, x0, lb, ub, rhobeg, rhoend, maxfun, npt, ftarget);
        else
            [x, fx, exitflag, nf, fhist, constrviolation, chist] = fbobyqa(fun, x0, lb, ub, rhobeg, rhoend, maxfun, npt, ftarget);
        end
    catch exception
        if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
            error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
        else % Private error; displayed as is
            rethrow(exception); 
        end
    end
end

% Postprocess the result 
try % postpdfo is a private function that may generate public errors; error-handeling needed
    [x, fx, exitflag, output] = postpdfo(x, fx, exitflag, output, nf, fhist, constrviolation, chist, options, probinfo);
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

% bobyqa ends
return
