function [x, fx, exitflag, output] = pdfo(varargin)
%PDFO is a package for solving the following generic continuous
%   optimization problem without using derivatives:
%
%   minimize    fun(x)
%       s.t.    Aineq * x <= bineq,
%               Aeq * x = beq,
%               lb <= x <= ub,
%               cineq(x) <= 0,
%               ceq(x) = 0.
%
%   In the backend, PDFO calls late Professor M.J.D. Powell's Fotran code 
%   UOBYQA, NEWUOA, BOBYQA, LINCOA, and COBYLA. 
%
%   The interface of PDFO is the same as that of function FMINCON included 
%   in the Optimization Toolbox of MATLAB. So PDFO can be called in the same 
%   way as calling FMINCON. In addition, PDFO can be called in some more 
%   flexible ways that are not allowed by FMINCON.
%
%   1. Basic syntax
%
%   The same as FMINCON, the command 
%    
%   x = pdfo(fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon) 
%
%   solves the problem formulated above, where
%
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
%   [x, fx, exitflag, output] = pdfo(INPUTS)
%
%   See "3. Outputs" below for explanations on these outputs. 
%
%   2. Flexible syntax
%
%   x = pdfo(fun, x0) solves 
%       minimize fun(x) 
%   x = pdfo(fun, x0, Aineq, bineq) solves
%       minimize fun(x) s.t. Aineq * x <= bineq 
%   x = pdfo(fun, x0, Aineq, bineq, Aeq, beq) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq
%   x = pdfo(fun, x0, Aineq, bineq, Aeq, beq, lb) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <= x
%   x = pdfo(fun, x0, Aineq, bineq, Aeq, beq, lb, ub) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <=x<= ub
%   x = pdfo(fun, x0, nonlcon) solves
%       minimize fun(x) s.t. cineq(x) <= 0, ceq(x) = 0
%   x = pdfo(fun, x0, Aineq, bineq, nonlcon) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, cineq(x) <= 0, ceq(x) = 0 
%   x = pdfo(fun, x0, Aineq, bineq, Aeq, beq, nonlcon) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, cineq(x) <= 0, ceq(x) = 0 
%   x = pdfo(fun, x0, Aineq, bineq, Aeq, beq, lb, nonlcon) solves
%       minimize fun(x) s.t. Aineq * x <= bineq, Aeq * x = beq, lb <= x, cineq(x) <= 0, ceq(x) = 0 
%
%   information = pdfo(request) returns information about the package
%       according to the information-requesting string "request", which can 
%       be 'about', 'author', 'email', 'maintainer', 'credits', 'copyright', 
%       'license', 'version', 'date', 'status', 'message', or 'information'.
%
%   3. Outputs
%
%   *** x is the approximate solution to the optimization problem
%   *** fx is fun(x)
%   *** exitflag is an integer indicating why PDFO or its backend solver 
%       returns; the possible values are  
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
%       code but cannot be returned by PDFO or its solvers
%   *** output is a structure with the following fields:
%       funcCount: number of function evaluations
%       nlcineq: cineq(x) (if there is nonlcon)
%       nlceq: ceq(x) (if there is nonlcon)
%       constrviolation: constrviolation of x (if problem is constrained)
%       fhist: history of function values
%       chist: history of constraint violations
%       solver: backend solver that does the computation
%       message: return message
%       warnings: a cell array that record all the warnings raised
%       during the computation
%   
%   4. Options 
%
%   The same as FMINCON, PDFO accepts options passed by a structure.
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
%       of the approximate solution obtained by PDFO; rhoend should be
%       positive and not larger than rhobeg; default: 1e-6
%   *** npt: number of interpolation points for constructing a model
%       (only for NEWUOA, BOBYQA, LINCOA); default value: 2*length(x0)+!
%   *** solver: a string indicating which solver to use; possible values are: 
%       'uobyqa', 'newuoa' (for unconstrained problems),
%       'bobyqa' (for bound-constrained or unconstrained problems),
%       'lincoa' (for linearly-constrained or bound-constrained or
%       unconstrained problems),
%       'cobyla' (for general constrained or unconstrained problems)
%   *** classical: a boolean value indicating whether to call the classical 
%       Powell code or not; default: false
%   *** scale: a boolean value that indicating whether to scale the problem
%       according to bounds or not; default: false
%   *** quiet: a boolean value indicating whether to keep quiet or not;
%       default: true (if false PDFO will print the return message of the
%       Fortran code)
%   *** debug: a boolean value indicating whether to debug or not; default: false
%   *** chkfunval: a boolean value indicating whether to verify the returned 
%       function and constraint (if applicable) value or not; default: false
%       (if true, PDFO will check whether the returned value of fun and nonlcon 
%       matches fun(x) and nonlcon(x) or not, which costs a function/constraint 
%       evaluation)
%
%   For example, the following code 
%
%   options = struct();
%   options.maxfun = 50;
%   x = pdfo(@cos, -1, 2, 3, options);
%
%   solves 
%       min cos(x) s.t. 2 * x <= 3 
%   starting from x0=-1 with at most 50 function evaluations.
%
%   5. Problem defined by a structure
%
%   The same as FMINCON, a problem can be passed to PDFO by a structure
%   PROBLEM containing the following fields: 
%   PROBLEM.objective, PROBLEM.x0, PROBLEM.Aineq, PROBLEM.bineq,
%   PROBLEM.Aeq, PROBLEM.beq, PROBLEM.lb, PROBLEM.ub, PROBLEM.nonlcon,
%   PROBLEM.options, where PROBLEM.objective is the function name or
%   function handle of the objective function (corresponding to the input
%   'fun' mentioned above), and all the other fields correspond to the
%   inputs introduced above with the same names. The backend solver can
%   be indicated by either PROBLEM.solver or PROBLEM.options.solver; if
%   both fields are defined, then PROBLEM.solver will be followed.
%
%   For example, the following code 
%
%   problem = struct();
%   problem.objective = @cos;
%   problem.x0 = -1;
%   problem.Aineq = 2;
%   problem.bineq = 3;
%   problem.options.maxfun = 50;
%   x = pdfo(problem);
%   
%   solves 
%       min cos(x) s.t. 2 * x <= 3 
%   starting from x0=-1 with at most 50 function evaluations.
%
%   See also UOBYQA, NEWUOA, BOBYQA, LINCOA, COBYLA.
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
% Remarks
%
% 1. Public function v.s. private function 
% 1.1. Public functions are functions that can be directly called by users.
% They should be either pdfo or a solver. 
% 1.2. Private functions are functions that are not supposed to be called
% by users. They are the preprocessing/postprocessing functions and
% auxiliary functions.
%
% 1. Errors that may be generated by pdfo
%
% 1.1. Normal error v.s. unexpected error
% A. Normal errors are usually caused by incorrect inputs.
% B. Unexpected errors usually imply bugs in the code. Such errors are 
% displayed with 'UNEXPECTED ERROR:', except those generated by mex files.
%
% 2.1. Public errors v.s. private errors
% A. Public errors are designed to be shown to users. The identifier
% of an public error starts with 'PUB_FUN:', where PUB_FUN is a 
% public function, i.e., the function that is being called by the user 
% (should be either pdfo or a solver). Thanks to the error-handling 
% code in public functions, public errors are displayed in a 'friendly' 
% manner as follows:
%
% Error using PUB_FUN (line XXX) 
% error message
% 
% where the error message is a string starts with 'PUB_FUN:'. If the
% error is not generated in the file of PUB_FUN, then an additional line
% as follows is displayed:
%
% (error generated in FILE, line YYY) 
%
% From the user's view point, a public error is an error raised by the
% function that he/she is calling. Therefore, if the user calls pdfo, 
% then PUB_FUN should be pdfo; if the user calls a solver directly, then 
% PUB_FUN should be the solver.
%
% B. Private errors are not expected to be seen by users unless there 
% is a bug. All private errors are unexpected errors (but not vice versa).
% The identifier of a private error starts with 'PRI_FUN:', where PRI_FUN 
% represents a private function. Private errors are displayed in the
% default manner, showing the detailed call stack trace, which is not so 
% friendly but easy to debug.
%
% C. Private functions display errors as is. Public functions catch 
% the exceptions thrown by private functions and display public/private
% errors properly.
%
% 2. Warnings that may be raised by pdfo
% 2.1. All the  warnings are considered to be public, i.e., they are
% desigened to be shown to the users.
% 2.2. Warnings are displayed without the call stack trace. 
%
% 3. probinfo 
% !!! TREAT probinfo AS A READONLY VARIABLE AFTER PREPDFO !!!
% !!! DO NOT CHANGE probinfo AFTER PREPDFO !!! 
%
% TODO:
% 1. Implicit NONE and variable declaration in the Fortran code
% 2. Change the interface of prepdfo to 
%    probinfo = prepdfo(argin{:}, interface_type),
%    where interface_type is one of 'unconstrained',
%    'bound-constrained', 'linearly-constrained',
%    'nonlinearly-constrained'.
%    All the information needed by the solvers should be included in probinfo.
% 3. To add a new solver, we only need to call prepdfo, call the solver
%    using the information in probinfo, record the results in output, and
%    then call postpdfo. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pdfo starts

callstack = dbstack;
funname = callstack(1).name; % Name of the current function 

% OUTPUT records the information that is produced by the solver and
% intended to pass to postpdfo.
% OUTPUT should contain at least x, fx, exitflag, funcCount, and constrviolation; 
% for internal solvers (solvers from PDFO), it should also contain fhist, chist, warnings;
% for lincoa, it should also contain constr_modified; 
% for nonlinearly constrained internal solvers, it should also contain nlcineq and nlceq. 
output = struct();

output.warnings = {}; % A cell that records all the warnings
% This version of pdfo.m produces no warning. However, initializing output.warnings 
% is still necessary, as output.warnings is required by postpdfo.
warning ('off', 'backtrace'); % Do not display the stack trace of a warning

maxarg = 10; % Maximal number of inputs
nvararg = length(varargin); % Number of inputs

% Interpret the input.
% Expected inputs: [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options],
% yet some of them may be omitted
if (nvararg < 1) % Public/normal error
    error(sprintf('%s:TooFewInputs', funname), '%s: at least 1 input.', funname);
elseif (nvararg == 1)
    % The only input should be either a information-requesting string or
    % a problem-defining structure.
    args = varargin; 
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
else % Public/normal error
    error(sprintf('%s:TooManyInputs', funname), '%s: at most %d inputs.', funname, maxarg);
end

% Preprocess the input
try % prepdfo and package_info are private functions that may generate public errors; error-handeling needed
    if (nvararg == 1) && (isa(varargin{1}, 'char') || isa(varargin{1}, 'string'))
        % If there is only 1 input and it is a string, then it should be
        % a string requesting information about the pacakge.
        x = package_info(varargin{1});
        return % Return immediately
    else
        [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options, probinfo] = prepdfo(args{:});
    end
catch exception 
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

if probinfo.infeasible % The problem turned out infeasible during prepdfo
    output.x = x0; 
    output.fx = fun(output.x);
    output.exitflag = -4;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_x0;
    output.chist = output.constrviolation;
    output.nlcineq = probinfo.nlcineq_x0;
    output.nlceq = probinfo.nlceq_x0;
    if strcmp(options.solver, 'lincoa') % LINCOA requires constr_modified to exist in output
        output.constr_modified = false;
    end
elseif probinfo.nofreex % x was fixed by the bound constraints during prepdfo
    output.x = probinfo.fixedx_value;
    output.fx = fun(output.x);
    output.exitflag = 13;
    output.funcCount = 1;
    output.fhist = output.fx;
    output.constrviolation = probinfo.constrv_fixedx;
    output.chist = output.constrviolation;
    output.nlcineq = probinfo.nlcineq_fixedx;
    output.nlceq = probinfo.nlceq_fixedx;
    if strcmp(options.solver, 'lincoa') % LINCOA requires constr_modified to exist in output
        output.constr_modified = false;
    end
else 
    % The problem turns out 'normal' during prepdfo. Solve it by
    % options.solver, which has been defined in prepdfo.
    try
        switch lower(options.solver)
        case 'uobyqa'
            [x, fx, exitflag, output] = uobyqa(fun, x0, options);
        case 'newuoa'
            [x, fx, exitflag, output] = newuoa(fun, x0, options);
        case 'bobyqa'
            [x, fx, exitflag, output] = bobyqa(fun, x0, lb, ub, options);
        case 'lincoa'
            [x, fx, exitflag, output] = lincoa(fun, x0, Aineq, bineq, Aeq, beq, lb, ub, options);
        case 'cobyla'
            [x, fx, exitflag, output] = cobyla(fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options);
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
end

% Postprocess the result 
try % postpdfo is a private function that may generate public errors; error-handeling needed
    [x, fx, exitflag, output] = postpdfo(probinfo, output);
catch exception
    if ~isempty(regexp(exception.identifier, sprintf('^%s:', funname), 'once')) % Public error; displayed friendly 
        error(exception.identifier, '%s\n(error generated in %s, line %d)', exception.message, exception.stack(1).file, exception.stack(1).line);
    else % Private error; displayed as is
        rethrow(exception); 
    end
end

% pdfo ends
return
