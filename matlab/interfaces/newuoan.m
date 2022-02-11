function [x, fx, exitflag, output] = newuoan(varargin)
%NEWUOA is a solver for solving the following unconstrained continuous
%   optimization problem without using derivatives:
%
%   minimize    fun(x).
%
%   In the backend, NEWUOA calls late Professor M.J.D. Powell's algorithm
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
%   *** exitflag is an integer indicating why NEWUOA returns; the possible values are
%       0: the lower bound for the trust region radius is reached
%       1: the target function value is achieved
%       2: a trust region step failed to reduce the quadratic model (possible only in classical mode)
%       3: the objective function has been evaluated maxfun times
%       14: a linear feasibility problem received and solved
%       20: the trust-region iteration has been performed for 10*maxfun times
%       -1: NaN occurs in x (possible only in the classical mode)
%       -2: the objective function returns an Inf/NaN value (possible only in classical mode)
%       -3: NaN occurs in the models (possible only in classical mode)
%   *** output is a structure with the following fields:
%       funcCount: number of function evaluations
%       xhist: history of iterates (if options.output_xhist = true)
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
%   *** fortran: a boolean value indicating whether to call Fortran code or
%       not; default: true
%   *** classical: a boolean value indicating whether to call the classical
%       version of Powell's Fortran code or not; default: false
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
%       -2, 3, or -3.
%       0: there will be no printing; this is the default;
%       1: a message will be printed to the screen at the return, showing
%          the best vector of variables found and its objective function value;
%       2: in addition to 1, at each "new stage" of the computation, a message
%          is printed to the screen with the best vector of variables so far
%          and its objective function value;
%       3: in addition to 2, each function evaluation with its variables will
%          be printed to the screen;
%       -1, -2, -3: the same information as 1, 2, 3 will be printed, not to
%          the screen but to a file named NEWUOA_output.txt; the file will be
%          created if it does not exist; the new output will be appended to
%          the end of this file if it already exists. Note that iprint = -3
%          can be costly in terms of time and space.
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
%       maxhist iterates of the algorithm; default: false
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
%   starting from x0 = -1 with at most 50 function evaluations.
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
%   starting from x0 = -1 with at most 50 function evaluations.
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
% Even if invoker = 'pdfon', we still need to call prepdfo, which will assign
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
npt = options.npt;
maxfun = options.maxfun;
rhobeg = options.rhobeg;
rhoend = options.rhoend;
eta1 = options.eta1;
eta2 = options.eta2;
gamma1 = options.gamma1;
gamma2 = options.gamma2;
ftarget = options.ftarget;
maxhist = options.maxhist;
output_xhist = options.output_xhist;
iprint = options.iprint;

if ~strcmp(invoker, 'pdfon') && probinfo.feasibility_problem
    % An "unconstrained feasibility problem" is rediculous, yet nothing wrong mathematically.
    output.x = x0;
    % We could set fx = [], funcCount = 0, and fhist = [] since no function evaluation
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
    % Call the Fortran code
    if options.classical
        fsolver = @fnewuoan_classical;
    else
        fsolver = @fnewuoan;
    end
    % The mexified Fortran Function is a private function generating only private errors;
    % however, public errors can occur due to, e.g., evalobj; error handling needed.
    try
        [x, fx, exitflag, nf, xhist, fhist] = ...
            fsolver(fun, x0, rhobeg, rhoend, eta1, eta2, gamma1, gamma2, ftarget, maxfun, npt, ...
            iprint, maxhist, double(output_xhist));
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
    if output_xhist
        output.xhist = xhist;
    end
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
