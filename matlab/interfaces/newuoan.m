function [x, fx, exitflag, output] = newuoan(varargin)
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
%   x = newuoan(fun, x0) 
%
%   solves the problem formulated above, where
%   *** fun is the name or function handle of the objective function
%   *** x0 is the starting point; x0 CANNOT be omitted or set to []
%
%   2. Outputs
%
%   The function can also be called with more outputs, e.g.,
%
%   [x, fx, exitflag, output] = newuoan(INPUTS)
%
%   *** x is the approximate solution to the optimization problem
%   *** fx is fun(x)
%   *** exitflag is an integer indicating why NEWUOA returns; the
%       possible values are 
%       0: the lower bound for the trust region radius is reached
%       1: the target function value is achieved
%       2: a trust region step failed to reduce the quadratic model
%       3: the objective function has been evaluated maxfun times
%       4, 7, 8, 9: rounding errors become severe in the Fortran code 
%       14: a linear feasibility problem received and solved
%       -1: NaN occurs in x
%       -2: the objective function returns an NaN or nearly infinite
%       value (only in the classical mode)
%       -3: NaN occurs in the models
%       exitflag = 5, 10, 11, 12 are possible exitflags of the Fortran
%       code but cannot be returned by NEWUOA 
%   *** output is a structure with the following fields:
%       funcCount: number of function evaluations
%       fhist: history of function values
%       solver: backend solver that does the computation, i.e., 'newuoan'
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
%   *** rhobeg: initial trust region radius; typically, rhobeg should be in
%       the order of one tenth of the greatest expected change to a variable;
%       rhobeg should be positive; default: 1 
%   *** rhoend: final trust region radius; rhoend reflects the precision
%       of the approximate solution obtained by NEWUOA; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** npt: number of interpolation points for constructing a model
%       default: 2*length(x0)+1
%   *** classical: a boolean value indicating whether to call the classical 
%       Powell code or not; default: false
%   *** quiet: a boolean value indicating whether to keep quiet or not;
%       default: true (if it is false, NEWUOA will print the return message of
%       the Fortran code)
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunval: a boolean value indicating whether to verify the returned 
%       function value or not; default: false
%       (if it is true, NEWUOA will check whether the returned value of fx
%       matches fun(x) or not, which costs a function evaluation;
%       designed only for debugging) 
%
%   For example, the following code 
%   
%   options = struct();
%   options.maxfun = 50;
%   x = newuoan(@cos, -1, options);
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
%   x = newuoan(problem);
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
% Remarks: 
% !!! TREAT probinfo AS A READONLY VARIABLE AFTER PREPDFO !!!
% !!! DO NOT CHANGE probinfo AFTER PREPDFO !!! 
%
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% newuoan starts
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
% Even if invoker='pdfon', we still need to call prepdfo, which will assign 
% values to fun, x0, ..., options.
try % prepdfo is a private function that may generate public errors; error-handling needed
    [fun, x0, ~, ~, ~, ~, ~, ~, ~, options, probinfo] = prepdfo(args{:}); 
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

% Extract the options
rhobeg = options.rhobeg;
rhoend = options.rhoend;
%eta1 = options.eta1;
%eta2 = options.eta2;
%gamma1 = options.gamma1;
%gamma2 = options.gamma2;
ftarget = options.ftarget;
maxfun = options.maxfun;
npt = options.npt;
%iprint = options.iprint;
% Here are the values of eta1, eta2, gamma1, and gamma2 used in Powell's
% code. They should be included in options later.
eta1 = 0.1;
eta2 = 0.7;
gamma1 = 0.5;
gamma2 = 2;
iprint = 0;

if ~strcmp(invoker, 'pdfon') && probinfo.feasibility_problem
    % An "unconstrained feasibility problem" is rediculous, yet nothing wrong mathematically.
    output.x = x0;
    % We could set fx=[], funcCount=0, and fhist=[] since no function evaluation 
    % occured. But then we will have to modify the validation of fx, funcCount, 
    % and fhist in postpdfo. To avoid such a modification, we set fx, funcCount, 
    % and fhist as below and then revise them in postpdfo.
    output.fx = fun(output.x);  % prepdfo has defined a fake objective function
    output.exitflag = 14;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = 0; % Unconstrained problem; set output.constrviolation to 0
    output.chist = []; % Unconstrained problem; set output.chist to []
else
    maxint = floor(0.99*min([gethuge('integer'), gethuge('mwSize'), gethuge('mwIndex')]));
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
    % The mexified Fortran Function is a private function generating only private errors; however, public errors can occur due to, e.g., evalobj; error handling needed 
        if options.classical
            [x, fx, exitflag, nf, fhist] = fnewuoan_classical(fun, x0, rhobeg, rhoend, maxfun, npt, ftarget);
        else
            [x, fx, exitflag, nf, fhist] = fnewuoan(fun, x0, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ftarget, maxfun, npt, iprint);
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
    output.constrviolation = 0; % Unconstrained problem; set output.constrviolation to 0
    output.chist = []; % Unconstrained problem; set output.chist to []
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

% newuoan ends
return
