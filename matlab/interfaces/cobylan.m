function [x, fx, exitflag, output] = cobylan(varargin)
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
%   In the backend, COBYLA calls late Professor M.J.D. Powell's algorithm
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
%   x = cobylan(fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon)
%
%   solves the problem formulated above, where
%   *** fun is the name or function handle of the objective function; if
%       there is no objective function (i.e., we have a feasibility problem),
%       then set fun = []
%   *** x0 is the starting point; x0 CANNOT be []
%   *** Aineq and bineq are the coefficient matrix and right-hand side of
%       the linear inequality constraint Aineq * x <= bineq; if there is
%       no such constraint, set Aineq = [], bineq = []
%   *** Aeq and beq are the coefficient matrix and right-hand side of the
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
%   [x, fx, exitflag, output] = cobylan(INPUTS)
%
%   See "3. Outputs" below for explanations on these outputs.
%
%   2. Flexible syntax
%
%   x = cobylan(fun, x0) solves
%       minimize fun(x)
%   x = cobylan(fun, x0, Aineq, bineq) solves
%       minimize fun(x) s.t. Aineq * x <= bineq
%   x = cobylan(fun, x0, Aineq, bineq, Aeq, beq) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq
%   x = cobylan(fun, x0, Aineq, bineq, Aeq, beq, lb) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <= x
%   x = cobylan(fun, x0, Aineq, bineq, Aeq, beq, lb, ub) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <= x <= ub
%   x = cobylan(fun, x0, nonlcon) solves
%       minimize fun(x) s.t. cineq(x) <= 0, ceq(x) = 0
%   x = cobylan(fun, x0, Aineq, bineq, nonlcon) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, cineq(x) <= 0, ceq(x) = 0
%   x = cobylan(fun, x0, Aineq, bineq, Aeq, beq, nonlcon) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, cineq(x) <= 0, ceq(x) = 0
%   x = cobylan(fun, x0, Aineq, bineq, Aeq, beq, lb, nonlcon) solves
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
%       14: a linear feasibility problem received and solved
%       15: a linear feasibility problem received but not solved
%       20: the trust-region iteration has been performed for 10*maxfun times
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
%       solver: backend solver that does the computation, i.e., 'cobylan'
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
%   *** ctol (only if classical = false; see below): tolerance for the constraint
%       validation; default: machine epsilon
%   *** rhobeg: initial trust region radius; typically, rhobeg should in the
%       order of one tenth of the greatest expected change to a variable;
%       rhobeg should be positive; default: 1 if the problem is not scaled,
%       0.5 if the problem is scaled
%   *** rhoend: final trust region radius; rhoend reflects the precision
%       of the approximate solution obtained by COBYLA; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** fortran: a boolean value indicating whether to call Fortran code or
%       not; default: true
%   *** classical: a boolean value indicating whether to call the classical
%       version of Powell's Fortran code or not; default: false
%   *** scale: a boolean value indicating whether to scale the problem
%       according to bounds or not; default: false; if the problem is to
%       be scaled, then rhobeg and rhoend mentioned above will be used as
%       the initial and final trust region radii for the scaled  problem
%   *** eta1, eta2, gamma1, gamma2 (only if classical = false)
%       eta1, eta2, gamma1, and gamma2 are parameters in the updating scheme
%       of the trust region radius. Roughly speaking, the trust region radius
%       is contracted by a factor of gamma1 when the reduction ratio is below
%       eta1, and  enlarged by a factor of gamma2 when the reduction ratio is
%       above eta2. It is required that 0 < eta1 <= eta2 < 1 and
%       0 < gamma1 < 1 < gamma2. Normally, eta1 <= 0.25. It is not recommended
%       to set eta1 >= 0.5. Default: eta1 = 0.1, eta2 = 0.7, gamma1 = 0.5,
%       and gamma2 = 2.
%   *** iprint: a flag deciding how much information will be printed during
%       the computation; possible values are value 0 (default), 1, -1, 2,
%       -2, 3, or -3:
%       0: there will be no printing;
%       1: a message will be printed to the screen at the return, showing
%          the best vector of variables found and its objective function value;
%       2: in addition to 1, at each "new stage" of the computation, a message
%          is printed to the screen with the best vector of variables so far
%          and its objective function value;
%       3: in addition to 2, each function evaluation with its variables will
%          be printed to the screen;
%       -1, -2, -3: the same information as 1, 2, 3 will be printed, not to
%          the screen but to a file named SOLVER_output.txt; the file will be
%          created if it does not exist; the new output will be appended to
%          the end of this file if it already exists. Note that iprint = -3
%          can be costly in terms of time and space.
%       When quiet = true (see below), setting iprint = 1, 2, or 3 is
%       the same as setting it to -1, -2, or -3, respectively.
%       Note:
%       When classical = true, only iprint = 0 is supported;
%       When fortran = true, only iprint = 0, -1, -2, -3 are supported
%       (due to I/O confliction between Fortran and MATLAB);
%       When quiet = true (see below), setting iprint = 1, 2, or 3 is
%       the same as setting it to -1, -2, or -3, respectively.
%   *** quiet: a boolean value indicating whether to keep quiet or not;
%       if this flag is set to false or not set, then it affects nothing;
%       if it is set to true and iprint = 1, 2, or 3, the effect is the
%       same as setting iprint to -1, -2, or -3, respectively; default: true
%   *** maxhist: a nonnegative integer controlling how much history will
%       be included in the output structure; default: maxfun;
%       *******************************************************************
%       IMPORTANT NOTICE:
%       If maxhist is so large that recording the history takes too much memory,
%       the Fortran code will reset maxhist to a smaller value. The maximal
%       amount of memory defined the Fortran code is 2GB.
%       *******************************************************************
%   *** output_xhist: a boolean value indicating whether to output the
%       history of the iterates; if it is set to true, then the output
%       structure will include a field "xhist", which contains the last
%       maxhist iterates of the algorithm; default: false;
%   *** output_nlchist: a boolean value indicating whether to output the
%       history of the function values; if it is set to true; then the
%       output structure will include fields "nlcihist" and "nlcehist",
%       which respectively contain the inequality and equality constraint
%       values of the last maxhist iterates of the algorithm; default: false
%   *** maxfilt: a nonnegative integer indicating the maximal length of the
%       "filter" used for selecting the returned solution; default: 2000
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunval: a boolean value indicating whether to verify the returned
%       function and constraint (if applicable) value or not; default: false
%       (if it is true, COBYLA will check whether the returned values of fun
%       and nonlcon matches fun(x) and nonlcon(x) or not, which costs
%       function/constraint evaluations; designed only for debugging)
%
%   For example, the following code
%
%   options = struct();
%   options.maxfun = 50;
%   x = cobylan(@cos, -1, 2, 3, options);
%
%   solves
%       min cos(x) s.t. 2 * x <= 3
%   starting from x0 = -1 with at most 50 function evaluations.
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
%   x = cobylan(problem);
%
%   solves
%       min cos(x) s.t. 2 * x <= 3
%   starting from x0 = -1 with at most 50 function evaluations.
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
% !!! TREAT probinfo AS A READONLY VARIABLE AFTER PREPDFO !!!
% !!! DO NOT CHANGE probinfo AFTER PREPDFO !!!
%
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cobylan starts

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
    % If 2 <= nvararg <= 10 and the last input is a structure or [], then it is the 'options'
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
% Even if invoker = 'pdfon', we still need to call prepdfo, which will assign
% values to fun, x0, ..., options.
try % prepdfo is a private function that may generate public errors; error-handling needed
    [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options, probinfo] = prepdfo(args{:});
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
    output.nlcineq = probinfo.nlcineq_x0;
    output.nlceq = probinfo.nlceq_x0;
elseif ~strcmp(invoker, 'pdfon') && probinfo.nofreex % x was fixed by the bound constraints during prepdfo
    output.x = probinfo.fixedx_value;
    output.fx = fun(output.x);
    output.exitflag = 13;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_fixedx;
    output.chist = output.constrviolation;
    output.nlcineq = probinfo.nlcineq_fixedx;
    output.nlceq = probinfo.nlceq_fixedx;
elseif ~strcmp(invoker, 'pdfon') && probinfo.feasibility_problem && ~strcmp(probinfo.refined_type, 'nonlinearly-constrained')
    output.x = x0;  % prepdfo has tried to set x0 to a feasible point (but may have failed)
    % We could set fx = [], funcCount = 0, and fhist = [] since no function evaluation
    % occured. But then we will have to modify the validation of fx, funcCount,
    % and fhist in postpdfo. To avoid such a modification, we set fx, funcCount,
    % and fhist as below and then revise them in postpdfo.
    output.fx = fun(output.x);  % prepdfo has defined a fake objective function
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_x0;
    output.chist = output.constrviolation;
    output.nlcineq = [];  % No nonlinear constraints
    output.nlceq = [];
    if output.constrviolation < eps  % Did prepdfo find a feasible point?
        output.exitflag = 14;
    else
        output.exitflag = 15;
    end
else % The problem turns out 'normal' during prepdfo
    % Include all the constraints into one single 'nonlinear constraint'
    funcon = @(x) cobylan_funcon(x, fun, Aineq, bineq, Aeq, beq, lb, ub, nonlcon);
    % Detect the number of the constraints (required by the Fortran code)
    [f_x0, constr_x0, m_nlcineq, m_nlceq] = funcon(x0);
    % m_nlcineq: number of nonlinear inequality constraints
    % m_nlceq: number of nonlinear equality constraints

    % If the constraint evaluation fails at x0, then m_nlcineq and m_nlceq are set to NaN by
    % funcon(x0). We have to raise an error and stop, because we do not know the number of
    % constraints, which is needed by the Fortran code.
    if (isnan(m_nlcineq) || isnan(m_nlceq))
        % Public/normal error
        error(sprintf('%s:ConstraintFailureAtX0', funname), '%s: The constraint evaluation fails at x0. %s cannot continue.', funname, funname);
    end

    if (length(constr_x0) > maxint())
        % Public/normal error
        error(sprintf('%s:ProblemTooLarge', funname), '%s: The problem is too large; at most %d constraints are allowed.', funname, maxint());
    end

    % Extract the options
    maxfun = options.maxfun;
    rhobeg = options.rhobeg;
    rhoend = options.rhoend;
    eta1 = options.eta1;
    eta2 = options.eta2;
    gamma1 = options.gamma1;
    gamma2 = options.gamma2;
    ftarget = options.ftarget;
    ctol = options.ctol;
    cweight = options.cweight;
    maxhist = options.maxhist;
    output_xhist = options.output_xhist;
    output_nlchist = options.output_nlchist;
    maxfilt = options.maxfilt;
    iprint = options.iprint;
    precision = options.precision;
    debug_flag = options.debug;
    if options.classical
        variant = 'classical';
    else
        variant = 'modern';
    end

    % Call the Fortran code
    solver = 'cobyla';
    fsolver = str2func(get_mexname(solver, precision, debug_flag, variant));
    % The mexified Fortran Function is a private function generating only private errors;
    % however, public errors can occur due to, e.g., evalobj; error handling needed.
    try
        [x, fx, constrviolation, constr, exitflag, nf, xhist, fhist, chist, conhist] = ...
            fsolver(funcon, x0, f_x0, constr_x0, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ...
            ftarget, ctol, cweight, maxfun, iprint, maxhist, double(output_xhist), ...
            double(output_nlchist), maxfilt);
        % Fortran MEX does not provide an API for reading Boolean variables. So we convert
        % output_xhist and output_nlchist to doubles (0 or 1) before passing them to the MEX gateway.
        % In C MEX, however, we have mxGetLogicals.
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
    if (output_xhist)
        output.xhist = xhist;
    end
    output.fhist = fhist;
    output.constrviolation = constrviolation;
    output.chist = chist;
    if output_nlchist
        output.nlcihist = -conhist(end-m_nlcineq-2*m_nlceq+1 : end-2*m_nlceq, :);
        if isempty(output.nlcihist)
            output.nlcihist = []; % We uniformly use [] to represent empty objects
        end
        output.nlcehist = -conhist(end-m_nlceq+1 : end, :);
        if isempty(output.nlcehist)
            output.nlcehist = []; % We uniformly use [] to represent empty objects
        end
    end
    % OUTPUT should also include nonlinear constraint values, if any
    output.nlcineq = [];
    output.nlceq = [];
    if ~isempty(nonlcon)
        output.nlcineq = -constr(end-m_nlcineq-2*m_nlceq+1 : end-2*m_nlceq);
        if isempty(output.nlcineq)
            output.nlcineq = []; % We uniformly use [] to represent empty objects
        end
        output.nlceq = -constr(end-m_nlceq+1 : end);
        if isempty(output.nlceq)
            output.nlceq = []; % We uniformly use [] to represent empty objects
        end
    end
end

% Postprocess the result
try % postdfo is a private function that may generate public errors; error-handling needed
    [x, fx, exitflag, output] = postpdfo(probinfo, output);
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception);
    end
end

% cobylan ends
return

%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%
function [f, constr, m_nlcineq, m_nlceq] = cobylan_funcon(x, fun, Aineq, bineq, Aeq, beq, lb, ub, nonlcon)
f = fun(x);
[constr, m_nlcineq, m_nlceq] = cobylan_con(x, Aineq, bineq, Aeq, beq, lb, ub, nonlcon);
return

function [constr, m_nlcineq, m_nlceq] = cobylan_con(x, Aineq, bineq, Aeq, beq, lb, ub, nonlcon)
% The Fortran backend takes at input a constraint: constr(x) >= 0
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
constr = [-cineq; ceq; -ceq];
if ~isempty(nonlcon)
    [nlcineq, nlceq, succ] = nonlcon(x); % Nonlinear constraints: nlcineq <= 0, nlceq = 0
    if succ
        m_nlcineq = length(nlcineq);
        m_nlceq = length(nlceq);
        constr = [constr; -nlcineq; nlceq; -nlceq];
    else
        % Evaluation of nonlcon fails.
        % In this case, we pass a SCALAR NaN to the MEX gateway, which will handle it properly.
        % Ideally, we should return an NaN vector with proper size, but the size is unknown here.
        m_nlcineq = NaN;
        m_nlceq = NaN;
        constr = NaN;
    end
else
    m_nlcineq = 0;
    m_nlceq = 0;
end
return
