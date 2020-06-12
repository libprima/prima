function [x, fx, exitflag, output] = bobyqan(varargin)
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
%   x = bobyqan(fun, x0, lb, ub) 
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
%   [x, fx, exitflag, output] = bobyqan(INPUTS)
%
%   See "3. Outputs" below for explanations on these outputs. 
%
%   2. Flexible syntax
%
%   x = bobyqan(fun, x0) solves 
%       minimize fun(x) 
%   x = bobyqan(fun, x0, lb) solves
%       minimize fun(x) s.t. lb <= x
%
%   3. Outputs
%
%   *** x is the approximate solution to the optimization problem
%   *** fx is fun(x)
%   *** exitflag is an integer indicating why BOBYQA returns; the
%       possible values are 
%       0: the lower bound for the trust region radius is reached
%       1: the target function value is achieved
%       2: a trust region step failed to reduce the quadratic model
%       3: the objective function has been evaluated maxfun times
%       4, 7, 8, 9: rounding errors become severe in the Fortran code 
%       13: all variables are fixed by the constraints
%       14: a linear feasibility problem received and solved
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
%       solver: backend solver that does the computation, i.e., 'bobyqan'
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
%   *** rhobeg: initial trust region radius; typically, rhobeg should be
%       in the order of one tenth of the greatest expected change to a variable;
%       rhobeg should be positive; default: min(1, min(ub-lb)/4) if the
%       problem is not scaled, 0.5 if the problem is scaled
%   *** rhoend: final trust region radius; rhoend reflects the precision
%       of the approximate solution obtained by BOBYQA; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** npt: number of interpolation points for constructing a model
%       default: 2*length(x0)+1
%   *** classical: a boolean value indicating whether to call the classical 
%       Powell code or not; default: false
%   *** scale: a boolean value indicating whether to scale the problem
%       according to bounds or not; default: false; if the problem is to be 
%       scaled, then rhobeg and rhoend mentioned above will be used as the 
%       initial and final trust region radii for the scaled  problem
%   *** honour_x0: a boolean value indicating whether to respect the
%       user-defiend x0 or not; default: false
%   *** quiet: a boolean value indicating whether to keep quiet or not;
%       default: true (if it is false, BOBYQA will print the return message of
%       the Fortran code)
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunval: a boolean value indicating whether to verify the returned 
%       function value or not; default: false
%       (if it is true, BOBYQA will check whether the returned value of fx
%       matches fun(x), which costs a function evaluation; designed only
%       for debugging) 
%
%   For example, the following code 
%   
%   options = struct();
%   options.maxfun = 50;
%   x = bobyqan(@cos, -1, 2, 3, options);
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
%   x = bobyqan(problem);
%   
%   solves 
%       min cos(x) s.t. 2 <= x <= 3 
%   starting from x0=-1 with at most 50 function evaluations.
%
%   See also PDFO, UOBYQA, NEWUOA, LINCOA, COBYLA.
%
%   See https://www.pdfo.net for more information.
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
% Remarks: 
% !!! TREAT probinfo AS A READONLY VARIABLE AFTER PREPDFO !!!
% !!! DO NOT CHANGE probinfo AFTER PREPDFO !!! 
%
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bobyqan starts

callstack = dbstack;
funname = callstack(1).name; % Name of the current function 
if length(callstack) >= 2
    invoker = callstack(2).name; % Name of the function who calls this function 
else
    invoker = '';
end
internal_invokers = {'pdfon'}; % Invokers from this package; may have others in the future

% OUTPUT records the information that is produced by the solver and
% intended to pass to postpdfo.
% OUTPUT should contain at least x, fx, exitflag, funcCount, and constrviolation; 
% for internal solvers (solvers from PDFO), it should also contain fhist, chist, warnings;
% for lincoa, it should also contain constr_modified; 
% for nonlinearly constrained internal solvers, it should also contain nlcineq and nlceq. 
output = struct();
% N.B.: DO NOT record anything in PROBINFO. If the solver is called by pdfo, 
% then postpdfo will do nothing; the real postprocessing will be done when 
% pdfo calls postpdfo using the OUTPUT returned by solver together with the 
% PROBINFO in pdfo; that said, in such a senario, the PROBINFO of this solver 
% will NOT be passed to the real postprocessing. Indeed, the PROBINFO of
% this solver is set to empty in prepdfo.

output.warnings = {}; % A cell that records all the warnings
warning('off', 'backtrace'); % Do not display the stack trace of a warning

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
% Even if invoker='pdfon', we still need to call prepdfo, which will assign 
% values to fun, x0, ..., options.
try % prepdfo is a private function that may generate public errors; error-handling needed
    [fun, x0, ~, ~, ~, ~, lb, ub, ~, options, probinfo] = prepdfo(args{:}); 
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

if ~strcmp(invoker, 'pdfon') && probinfo.infeasible % The problem turned out infeasible during prepdfo
    output.x = x0; 
    output.fx = fun(output.x);
    output.exitflag = -4;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_x0;
    output.chist = output.constrviolation;
elseif ~strcmp(invoker, 'pdfon') && probinfo.nofreex % x was fixed by the bound constraints during prepdfo
    output.x = probinfo.fixedx_value;
    output.fx = fun(output.x);
    output.exitflag = 13;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_fixedx;
    output.chist = output.constrviolation;
elseif ~strcmp(invoker, 'pdfon') && probinfo.feasibility_problem
    % A "feasibility problem" with only bound constraints is rediculous yet nothing wrong mathematically
    output.x = x0;  % prepdfo has set x0 to a feasible point
    % We could set fx=[], funcCount=0, and fhist=[] since no function evaluation 
    % occured. But then we will have to modify the validation of fx, funcCount, 
    % and fhist in postpdfo. To avoid such a modification, we set fx, funcCount, 
    % and fhist as below and then revise them in postpdfo.
    output.fx = fun(output.x);  % prepdfo has defined a fake objective function
    output.exitflag = 14;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_x0;
    output.chist = output.constrviolation;
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
        if strcmp(invoker, 'pdfon')
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
    try % The mexified Fortran function is a private function generating only private errors; however, public errors can occur due to, e.g., evalobj; error handling needed 
        if options.classical 
            [x, fx, exitflag, nf, fhist, constrviolation, chist] = fbobyqan_classical(fun, x0, lb, ub, rhobeg, rhoend, maxfun, npt, ftarget);
        else
            [x, fx, exitflag, nf, fhist, constrviolation, chist] = fbobyqan(fun, x0, lb, ub, rhobeg, rhoend, maxfun, npt, ftarget);
        end
    catch exception
        if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
            error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
        else % Private error; displayed as is
            rethrow(exception); 
        end
    end
    % Record the results of the solver in OUTPUT
    output.x = x;
    output.fx = fx;
    output.exitflag = exitflag;
    output.funcCount = nf;
    output.fhist = fhist;
    output.constrviolation = constrviolation;
    output.chist = chist;
end

% Postprocess the result 
try % postpdfo is a private function that may generate public errors; error-handling needed
    [x, fx, exitflag, output] = postpdfo(probinfo, output);
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

% bobyqan ends
return
