function [x, fx, exitflag, output] = cobyla(varargin)
%COBYLA is a solver for solving the following generic continuous
%   optimization problem without using derivatives:
%
%   minimize    fun(x)
%       s.t.    Aineq * x <= bineq,
%               Aeq * x = beq,
%               lb <= x <= ub,
%               cineq(x) <= 0,
%               ceq(x) = 0.
%
%   In the backend, COBYLA calls late Professor M.J.D. Powell's Fotran code 
%   with the same name. The algorithm is described in [M. J. D. Powell,
%   A direct search optimization method that models the objective and 
%   constraint functions by linear interpolation, In Advances in Optimization 
%   and Numerical Analysis, eds. S. Gomez and J. P. Hennart, pages 51--67, 
%   Springer Verlag, Dordrecht, Netherlands, 1994].
%
%   The interface of COBYLA is the same as that of function FMINCON included 
%   in the Optimization Toolbox of MATLAB. So COBYLA can be called in the same 
%   way as calling FMINCON. In addition, COBYLA can be called in some more 
%   flexible ways that are not allowed by FMINCON.
%
%   1. Basic syntax
%
%   The same as FMINCON, the command
%
%   x = cobyla(fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon) 
%
%   solves the problem formulated above, where
%   *** fun is the name or function handle of the objective function; if
%       there is no objective function (i.e., we have a feasibility problem), 
%       then set fun = []
%   *** x0 is the starting point; x0 CANNOT be []
%   *** Aineq and bineq are the coeffcient matrix and right-hand side of 
%       the linear inequality constraint Aineq * x <= bineq; if there is
%       no such constraint, set Aineq = [], bineq = []
%   *** Aeq and beq are the coeffcient matrix and right-hand side of the 
%       linear equality constraint Aeq * x = beq; if there is no such
%       constraint, set Aeq = [], beq = []
%   *** lb and ub, which are vectors of the same length as x, are the 
%       lower and upper bound in the bound constraint lb <= x <= ub; 
%       set lb = [] if no lower bound, and ub = [] if no upper bound 
%   *** nonlcon is a function that has 1 input x and 2 outputs [cineq, ceq]; 
%       it calculates cineq(x) and ceq(x) for any given x; if the first
%       output of nonlcon is [], then there is no inequality constraint
%       cineq(x) <= 0; if the second output of nonlcon is [], then there
%       is no equality constraint ceq(x) = 0. If there is no nonlinear
%       constraint, set nonlcon = []
%
%   The function can also be called with more outputs, e.g.,
%
%   [x, fx, exitflag, output] = cobyla(INPUTS)
%
%   See "3. Outputs" below for explanations on these outputs. 
%
%   2. Flexible syntax
%
%   x = cobyla(fun, x0) solves 
%       minimize fun(x) 
%   x = cobyla(fun, x0, Aineq, bineq) solves
%       minimize fun(x) s.t. Aineq * x <= bineq 
%   x = cobyla(fun, x0, Aineq, bineq, Aeq, beq) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq
%   x = cobyla(fun, x0, Aineq, bineq, Aeq, beq, lb) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <= x
%   x = cobyla(fun, x0, Aineq, bineq, Aeq, beq, lb, ub) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <=x<= ub
%   x = cobyla(fun, x0, nonlcon) solves
%       minimize fun(x) s.t. cineq(x) <= 0, ceq(x) = 0
%   x = cobyla(fun, x0, Aineq, bineq, nonlcon) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, cineq(x) <= 0, ceq(x) = 0 
%   x = cobyla(fun, x0, Aineq, bineq, Aeq, beq, nonlcon) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, cineq(x) <= 0, ceq(x) = 0 
%   x = cobyla(fun, x0, Aineq, bineq, Aeq, beq, lb, nonlcon) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <= x, cineq(x) <= 0, ceq(x) = 0 
%
%   3. Outputs
%
%   *** x is the approximate solution to the optimization problem
%   *** fx is fun(x)
%   *** exitflag is an integer indicating why COBYLA returns; the
%       possible values are 
%       0: the lower bound for the trust region radius is reached
%       1: the target function value is achieved
%       2: a trust region step failed to reduce the quadratic model
%       3: the objective function has been evaluated maxfun times
%       4, 7, 8, 9: rounding errors become severe in the Fortran code 
%       13: all variables are fixed by the constraints
%       -1: NaN occurs in x
%       -2: the objective/constraint function returns NaN or nearly
%       infinite values (only in the classical mode)
%       -3: NaN occurs in the models
%       -4: constraints are infeasible 
%       exitflag = 5, 10, 11, 12 are possible exitflags of the Fortran
%       code but cannot be returned by COBYLA 
%   *** output is a structure with the following fields:
%       funcCount: number of function evaluations
%       nlcineq: cineq(x) (if there is nonlcon)
%       nlceq: ceq(x) (if there is nonlcon)
%       constrviolation: constrviolation of x (if problem is constrained)
%       fhist: history of function values
%       chist: history of constraint violations
%       solver: backend solver that does the computation, i.e., 'cobyla'
%       message: return message
%       warnings: a cell array that records all the warnings raised
%       during the computation
%   
%   4. Options 
%
%   The same as FMINCON, COBYLA accepts options passed by a structure.
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
%       of the approximate solution obtained by COBYLA; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** classical: a boolean value indicating whether to call the classical 
%       Powell code or not; default: false
%   *** scale: a boolean value that indicating whether to scale the problem
%       according to bounds or not; default: false
%   *** quiet: a boolean value indicating whether to keep quiet or not;
%       default: true (if false COBYLA will print the return message of the
%       Fortran code)
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunval: a boolean value indicating whether to verify the returned 
%       function and constraint (if applicable) value or not; default: false
%       (if true, COBYLA will check whether the returned value of fun and nonlcon 
%       matches fun(x) and nonlcon(x) or not, which costs a function/constraint 
%       evaluation) 
%
%   For example, the following code 
%
%   options = struct();
%   options.maxfun = 50;
%   x = cobyla(@cos, -1, 2, 3, options);
%
%   solves 
%       min cos(x) s.t. 2 * x <= 3 
%   starting from x0=2 with at most 50 function evaluations.
%
%   5. Problem defined by a structure
%
%   The same as FMINCON, a problem can be passed to COBYLA by a structure
%   PROBLEM containing the following fields: 
%   PROBLEM.objective, PROBLEM.x0, PROBLEM.Aineq, PROBLEM.bineq,
%   PROBLEM.Aeq, PROBLEM.beq, PROBLEM.lb, PROBLEM.ub, PROBLEM.nonlcon,
%   PROBLEM.options, where PROBLEM.objective is the function name or
%   function handle of the objective function (corresponding to the input 
%   'fun' mentioned above), and all the other fields correspond to the
%   inputs introduced above with the same names. 
%
%   For example, the following code 
%
%   problem = struct();
%   problem.objective = @cos;
%   problem.x0 = -1;
%   problem.Aineq = 2;
%   problem.bineq = 3;
%   problem.options.maxfun = 50;
%   x = cobyla(problem);
%   
%   solves 
%       min cos(x) s.t. 2 * x <= 3 
%   starting from x0=-1 with at most 50 function evaluations.
%
%   See also PDFO, UOBYQA, NEWUOA, BOBYQA, LINCOA.
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
% !!! TREAT probinfo and options AS READONLY VARIABLES AFTER PREPDFO!!!
% !!! DO NOT MODIFY THE INFORMATION IN probinfo OR options AFTER PREPDFO !!! 
%
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cobyla starts 

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

maxarg = 10; % Maximal number of inputs
nvararg = length(varargin); % Number of inputs

% Interpret the input.
% Expected inputs: [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options],
% yet some of them may be ommitted
if (nvararg < 1)
    if ismember(invoker, internal_invokers) % Private/unexpected error
        error(sprintf('%s:TooFewInputs', funname), '%s: UNEXPECTED ERROR: at least 1 input.', funname);
    else % Public/normal error
        error(sprintf('%s:TooFewInputs', funname), '%s: at least 1 input.', funname);
    end
elseif (nvararg == 1)
    args = varargin; % If there is only 1 input, then it is a structure specifying the problem
elseif (nvararg >= 2 && nvararg <= maxarg) 
    % If 2<=nvararg<=10 and the last input is a structure or [], then it is the 'options'
    if isempty(varargin{end}) || isa(varargin{end}, 'struct') 
        % If nvararg >= 4 and the second last input is a function, then it is the 'nonlcon'
        if (nvararg >= 4) && (isa(varargin{end-1}, 'char') || isa(varargin{end-1}, 'string') || isa(varargin{end-1}, 'function_handle'))
            args = [varargin(1:end-2), cell(1, maxarg-nvararg), varargin(end-1:end)]; % 'augment' the inputs to 10 by adding []; args{:} (should have 10 entries) will be the inputs for prepdfo
            % cell(m,n) returns an mxn array of [] 
        else
            args = [varargin(1:end-1), cell(1, maxarg-nvararg), varargin(end)];
        end
    % if nvararg >= 3 and the last input is a function, then it is the 'nonlcon'
    elseif (nvararg >= 3) && (isa(varargin{end}, 'char') || isa(varargin{end}, 'string') || isa(varargin{end}, 'function_handle'))
        args = [varargin(1:end-1), cell(1, maxarg-nvararg-1), varargin(end), {[]}];
    else
        args = [varargin, cell(1, maxarg-nvararg)];
    end
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
    [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options, probinfo] = prepdfo(args{:});
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

if ~strcmp(invoker, 'pdfo') && probinfo.infeasible % The problem turned out infeasible during prepdfo
    output.x = x0; 
    output.fx = fun(output.x);
    output.exitflag = -4;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_x0;
    output.chist = output.constrviolation;
    output.nlcineq = probinfo.nlcineq_x0;
    output.nlceq = probinfo.nlceq_x0;
elseif ~strcmp(invoker, 'pdfo') && probinfo.nofreex % x was fixed by the bound constraints during prepdfo
    output.x = probinfo.fixedx_value;
    output.fx = fun(output.x);
    output.exitflag = 13;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_fixedx;
    output.chist = output.constrviolation;
    output.nlcineq = probinfo.nlcineq_fixedx;
    output.nlceq = probinfo.nlceq_fixedx;
else % The problem turns out 'normal' during prepdfo
    % Include all the constraints into one single 'nonlinear constraint'
    con = @(x) cobyla_con(x, Aineq, bineq, Aeq, beq, lb, ub, nonlcon);
    % Detect the number of the constraints (required by the Fortran code)
    [conval_x0, m_nlcineq, m_nlceq] = con(x0);
    % m_nlcineq = number of nonlinear inequality constraints
    % m_nlceq = number of nonlinear equality constraints
    m = length(conval_x0); % number of constraints
    % Extract the options
    maxfun = options.maxfun;
    rhobeg = options.rhobeg;
    rhoend = options.rhoend;
    ftarget = options.ftarget;
    
    % Check whether the problem is too large for the Fortran code.
    % In the mex gateway, a workspace of size 
    % nw = n*(3*n+2*m+11)+4*m+6
    % will be allocated, which is the largest memory allocated by
    % COBYLA. If the value assigned to nw is so large that overflow
    % occurs, then there will be a Segmentation Fault!!! 
    % The largest possible value of nw depends on the type of nw in the
    % mex file, which is the default INTEGER type (~2E9 for integer*4, 
    % and ~9E18 for integer*8). This imposes an upper limit on the size
    % of problem solvable by this code. If nw is INTEGER*4, assuming
    % that m=10n, the largest value of n is ~9600. COBYLA is not
    % designed for so large problems. 
    % In the following code, gethuge returns the largest possible value
    % of the given data type in the mex environment.
    
    % The largest integer in the mex functions; the factor 0.99 provides a buffer
    maxint = floor(0.99*min([gethuge('integer'), gethuge('mwSize'), gethuge('mwIndex')]));
    n = length(x0);
    nw = n*(3*n+2*m+11)+4*m+6;
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
        wid = sprintf('%s:MaxfunTooLarge', funname);
        wmessage = sprintf('%s: maxfun exceeds the upper limit of Fortran integers; it is set to %d.', funname, maxfun);
        warning(wid, '%s', wmessage);
        output.warnings = [output.warnings, wmessage];
    end

    % Call the Fortran code
    try % The mexified Fortran function is a private function generating only private errors; however, public errors can occur due to, e.g., evalobj and evalcon; error handeling needed 
        if options.classical
            [x, fx, exitflag, nf, fhist, conval, constrviolation, chist] = fcobyla_classical(fun, con, x0, rhobeg, rhoend, maxfun, m, ftarget, conval_x0);
        else
            [x, fx, exitflag, nf, fhist, conval, constrviolation, chist] = fcobyla(fun, con, x0, rhobeg, rhoend, maxfun, m, ftarget, conval_x0);
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
    % OUTPUT should also include nonlinear constraint values, if any
    output.nlcineq = [];
    output.nlceq = [];
    if ~isempty(nonlcon)
        output.nlcineq = -conval(end-m_nlcineq-2*m_nlceq+1 : end-2*m_nlceq);
        if isempty(output.nlcineq) 
            output.nlcineq = []; % We uniformly use [] to represent empty objects
        end
        output.nlceq = -conval(end-m_nlceq+1 : end);
        if isempty(output.nlceq)
            output.nlceq = []; % We uniformly use [] to represent empty objects
        end
    end
end

% Postprocess the result 
try % postdfo is a private function that may generate public errors; error-handeling needed
    [x, fx, exitflag, output] = postpdfo(probinfo, options, output);
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

% cobyla ends
return

%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%
function [conval, m_nlcineq, m_nlceq] = cobyla_con(x, Aineq, bineq, Aeq, beq, lb, ub, nonlcon) 
% The Fortran backend takes at input a constraint: conval(x) >= 0 
% m_nlcineq = number of nonlinear inequality constraints
% m_nlceq = number of nonlinear equality constraints
cineq = [lb(lb>-inf) - x(lb>-inf); x(ub<inf) - ub(ub<inf)]; 
ceq = [];
if ~isempty(Aineq)
    cineq = [Aineq*x - bineq; cineq];
end
if ~isempty(Aeq)
    ceq = [ceq; Aeq*x - beq];
end
conval = [-cineq; ceq; -ceq];
if ~isempty(nonlcon)
    [nlcineq, nlceq] = nonlcon(x); % Nonlinear constraints: nlcineq <= 0, nlceq = 0
    m_nlcineq = length(nlcineq);
    m_nlceq = length(nlceq);
    conval = [conval; -nlcineq; nlceq; -nlceq];
else
    m_nlcineq = 0;
    m_nlceq = 0;
end
return
