function [x, fx, exitflag, output] = lincoan(varargin)
%LINCOA is a solver for solving the following linearly-constrained  continuous
%   optimization problem without using derivatives:
%
%   minimize    fun(x)
%       s.t.    Aineq * x <= bineq,
%               Aeq * x = beq,
%               lb <= x <= ub.
%
%   In the backend, LINCOA calls late Professor M.J.D. Powell's algorithm
%   with the same name.
%
%   The interface of LINCOA is similar to that of function FMINCON included
%   in the Optimization Toolbox of MATLAB, except that LINCOA does not accept
%   nonlinear constraints. So LINCOA can be called in a way similar to calling
%   FMINCON. In addition, LINCOA can also be called in some more flexible ways.
%
%   1. Basic syntax
%
%   Similar to FMINCON, the command
%
%   x = lincoan(fun, x0, Aineq, bineq, Aeq, beq, lb, ub)
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
%
%   The function can also be called with more outputs, e.g.,
%
%   [x, fx, exitflag, output] = lincoan(INPUTS)
%
%   See "3. Outputs" below for explanations on these outputs.
%
%   2. Flexible syntax
%
%   x = lincoan(fun, x0) solves
%       minimize fun(x)
%   x = lincoan(fun, x0, Aineq, bineq) solves
%       minimize fun(x) s.t. Aineq * x <= bineq
%   x = lincoan(fun, x0, Aineq, bineq, Aeq, beq) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq
%   x = lincoan(fun, x0, Aineq, bineq, Aeq, beq, lb) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <= x
%
%   3. Outputs
%
%   *** x is the approximate solution to the optimization problem
%   *** fx is fun(x)
%   *** exitflag is an integer indicating why LINCOA returns; the
%       possible values are
%       0: the lower bound for the trust region radius is reached
%       1: the target function value is achieved
%       2: a trust region step failed to reduce the quadratic model
%       3: the objective function has been evaluated maxfun times
%       4, 7, 8, 9: rounding errors become severe in the Fortran code
%       13: all variables are fixed by the constraints
%       14: a feasibility problem received and solved
%       15: a feasibility problem received but not solved
%       20: the trust-region iteration has been performed for 10*maxfun times
%       -1: NaN occurs in x
%       -2: the objective function returns an NaN or nearly infinite
%       value (only in the classical mode)
%       -3: NaN occurs in the models
%       -4: constraints are infeasible
%       exitflag = 5, 10, 11, 12 are possible exitflags of the Fortran
%       code but cannot be returned by LINCOA
%   *** output is a structure with the following fields:
%       funcCount: number of function evaluations
%       constrviolation: constrviolation of x (if problem is constrained)
%       fhist: history of function values
%       chist: history of constraint violations
%       solver: backend solver that does the computation, i.e., 'lincoan'
%       message: return message
%       warnings: a cell array that records all the warnings raised
%       during the computation
%
%   4. Options
%
%   The same as FMINCON, LINCOA accepts options passed by a structure.
%   Such a structure should be passed as an additional input appended to
%   the end of the input list in the basic syntax or the flexible syntax.
%
%   The options include
%   *** maxfun: maximal number of function evaluations; default: 500*length(x0)
%   *** ftarget: target function value; default: -Inf
%   *** ctol (only if classical = false; see below): tolerance for the constraint
%       validation; default: machine epsilon
%   *** rhobeg: initial trust region radius; typically, rhobeg should be in
%       the order of one tenth of the greatest expected change to a variable;
%       rhobeg should be positive; default: 1 if the problem is not scaled,
%       0.5 if the problem is scaled
%   *** rhoend: final trust region radius; rhoend reflects the precision
%       of the approximate solution obtained by LINCOA; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** npt: number of interpolation points for constructing a model
%       default: 2*length(x0)+1
%   *** fortran: a boolean value indicating whether to call Fortran code or
%       not; default: true
%   *** classical: a boolean value indicating whether to call the classical
%       version of Powell's Fortran code or not; default: false
%   *** scale: a boolean value indicating whether to scale the problem
%       according to bounds or not; default: false; if the problem is to be
%       scaled, then rhobeg and rhoend mentioned above will be used as the
%       initial and final trust region radii for the scaled  problem
%   *** iprint: a flag deciding how much information will be printed during
%       the computation; possible values are value 0 (default), 1, -1, 2,
%       -2, 3, or -3:
%       0: there will be no printing;
%       1: a message will be printed to the screen at the return, showing
%          the best vector of veriables found and its objective function value;
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
%   *** maxfilt: a nonnegative integer indicating the maximal length of the
%       "filter" used for selecting the returned solution; default: 2000
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunvsl: a boolean value indicating whether to verify the returned
%       function value or not; default: false
%       (if it is true, LINCOA will check whether the returned value of fun
%       matches fun(x) or not, which costs a function evaluation;
%       designed only for debugging)
%
%   For example, the following code
%
%   options = struct();
%   options.maxfun = 50;
%   x = lincoan(@cos, -1, 2, 3, options);
%
%   solves
%       min cos(x) s.t. 2 * x <= 3
%   starting from x0 = -1 with at most 50 function evaluations.
%
%   5. Problem defined by a structure
%
%   The same as FMINCON, a problem can be passed to LINCOA by a structure
%   PROBLEM containing the following fields:
%   PROBLEM.objective, PROBLEM.x0, PROBLEM.Aineq, PROBLEM.bineq,
%   PROBLEM.Aeq, PROBLEM.beq, PROBLEM.lb, PROBLEM.ub, POBLEM.options,
%   where PROBLEM.objective is the function name or function handle of
%   the objective function (corresponding to the input 'fun' mentioned
%   above), and all the other fields correspond to the inputs introduced
%   above with the same names.
%
%   For example, the following code
%
%   problem = struct();
%   problem.objective = @cos;
%   problem.x0 = -1;
%   problem.Aineq = 2;
%   problem.bineq = 3;
%   problem.options.maxfun = 50;
%   x = lincoan(problem);
%
%   solves
%       min cos(x) s.t. 2 * x <= 3
%   starting from x0 = -1 with at most 50 function evaluations.
%
%   See also PDFO, UOBYQA, NEWUOA, BOBYQA, COBYLA.
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

% lincoan starts

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
% for lincoan, it should also contain constr_modified;
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

maxarg = 9; % Maximal number of inputs
nvararg = length(varargin); % Number of inputs

% Interpret the input.
% Expected inputs: [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, options],
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
    % If 2 <= nvararg <= 9 and the last input is a structure (or []), then it is the 'options'
    if isa(varargin{end}, 'struct')
        varargin = [varargin(1:end-1), cell(1, maxarg-nvararg), varargin(end)]; % 'augment' the inputs to maxarg by adding []
        % cell(m,n) returns an mxn array of []
    else
        varargin = [varargin, cell(1, maxarg-nvararg)];  % 'augment' the inputs to maxarg by adding []
    end
    args = [varargin(1:8), {[]}, varargin(end)]; % args{:} (should have 10 entries) will be the inputs for prepdfo
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
    [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, ~, options, probinfo] = prepdfo(args{:});
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
    output.constr_modified = false;
elseif ~strcmp(invoker, 'pdfon') && probinfo.nofreex % x was fixed by the bound constraints during prepdfo
    output.x = probinfo.fixedx_value;
    output.fx = fun(output.x);
    output.exitflag = 13;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_fixedx;
    output.chist = output.constrviolation;
    output.constr_modified = false;
elseif ~strcmp(invoker, 'pdfon') &&  probinfo.feasibility_problem
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
    output.constr_modified = false; % LINCOA requires constr_modified to exist in output
    if output.constrviolation < eps  % Did prepdfo find a feasible point?
        output.exitflag = 14;
    else
        output.exitflag = 15;
    end
else % The problem turns out 'normal' during prepdfo
    % Include all the constraints into one single linear constraint
    % (A_aug)'*x <= b_aug; note the TRANSPOSE due to the data structure of
    % the Fortran code.
    n = length(x0);
    idmatrix = eye(n, n); % If we want to use sparse matrix, the mex gateway has to be modified
    Alb = idmatrix(lb > -inf, :);
    Aub = idmatrix(ub < inf, :);
    lb = lb(lb > -inf); % Remove infinity bounds!
    ub = ub(ub < inf); % Remove infinity bounds!
    % Note that prepdfo has removed the 'trivial' linear constraints in
    % [Aineq, bineq] and [Aeq, beq].
    A_aug = [Aineq; Aeq; -Aeq; -Alb; Aub]';
    b_aug = [bineq(:); beq(:); -beq(:); -lb(:); ub(:)];
    if ~(isempty(A_aug) && isempty(b_aug)) && ~(size(A_aug,1) == n && size(A_aug,2) == length(b_aug))
        % Private/unexpected error
        error(sprintf('%s:InvalidAugLinCon', funname), '%s: UNEXPECTED ERROR: invalid augmented linear constraints.', funname);
    end
    if isempty(A_aug)
    % We uniformly use [] to represent empty objects; its size is 0x0.
    % Changing this may cause matrix dimension inconsistency
        A_aug = [];
        b_aug = [];
    end

    if (length(b_aug) > maxint())
        % Public/normal error
        error(sprintf('%s:ProblemTooLarge', funname), '%s: The problem is too large; at most %d constraints are allowed.', funname, maxint());
    end

    % Extract the options
    npt = options.npt;
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
    maxfilt = options.maxfilt;
    iprint = options.iprint;
    precision = options.precision;
    debug_flag = options.debug;
    if options.classical
        variant = 'classical';
    else
        variant = 'modern';
    end

    % If x0 is not feasible, LINCOA will modify the constraints to make
    % it feasible (which is a bit strange).
    % prepdfo has tried to make find a feasible x0. Raise a warning is
    % x0 is not 'feasible enough' so that the constraints will be modified.
    if ~isempty(A_aug) && any(A_aug'*x0 > b_aug + 1e-10*max(1, abs(b_aug)))
        output.constr_modified = true;
        wid = sprintf('%s:InfeasibleX0', funname);
        wmsg = sprintf('%s: preprocessing code did not find a feasible x0; problem is likely infeasible; %s will modify the right-hand side of the constraints to make x0 feasible.', funname, funname);
        warning(wid, '%s', wmsg);
        output.warnings = [output.warnings, wmsg];
    else
        output.constr_modified = false;
    end

    % Call the Fortran code
    solver = 'lincoa';
    fsolver = str2func(get_mexname(solver, precision, debug_flag, variant));
    % The mexified Fortran function is a private function generating only private errors;
    % however, public errors can occur due to, e.g., evalobj; error handling needed
    try
        [x, fx, constrviolation, exitflag, nf, xhist, fhist, chist] = ...
            fsolver(fun, x0, A_aug, b_aug, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ...
            ftarget, ctol, cweight, maxfun, npt, iprint, maxhist, double(output_xhist), maxfilt);
        % Fortran MEX does not provide an API for reading Boolean variables. So we convert
        % output_xhist to a double (0 or 1) before passing it to the MEX gateway.
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

% lincoan ends
return
