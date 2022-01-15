function [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options, probinfo] = prepdfo(fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options)
%PREPDFO preprocesses the input to pdfo and its solvers.
%
%   ***********************************************************************
%   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%
%   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: private (not supposed to be called by users)
%
% Remarks
% 1. Input/output names: MATLAB allows to use the same name for inputs and outputs.
% 2. invoker: invoker is the function that calls prepdfo
%
% TODO: None
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepdfo starts

warnings = {}; % A cell that records all the warnings, will be recorded in probinfo

% Who is calling this function? Is it a correct invoker?
invoker_list = {'pdfon', 'uobyqan', 'newuoan', 'bobyqan', 'lincoan', 'cobylan'};
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if (length(callstack) == 1 || ~ismember(callstack(2).name, invoker_list))
    % Private/unexpected error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s should only be called by %s.', funname, funname, mystrjoin(invoker_list, ', '));
else
    invoker = callstack(2).name; % Name of the function who calls this function
end

if (nargin ~= 1) && (nargin ~= 10)
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: 1 or 10 inputs.', funname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If invoker is a solver called by pdfo, then prepdfo should have been called in pdfo.
if (length(callstack) >= 3) && strcmp(callstack(3).name, 'pdfon')
    if nargin ~= 10 % There should be 10 input arguments
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), ...
        '%s: UNEXPECTED ERROR: %d inputs received; this should not happen as prepdfo has been called once in pdfo.', funname, nargin);
    end
    % In this case, we set probinfo to empty.
    probinfo = [];
    return % Return because prepdfo has already been called in pdfo.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decode the problem if it is defined by a structure.
if (nargin == 1)
    [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options, warnings] = decode_problem(invoker, fun, warnings);
end

% Save problem information in probinfo.
% At return, probinfo has the following fields:
% 1. raw_data: problem data before preprocessing/validating, including
%    fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options.
%    raw_data is set to struct() unless in debug mode.
% 2. refined_data: problem data after preprocessing/validating, including
%    fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options.
%    refined_data is set to struct() unless in debug mode or the problem is scaled.
% 3. fixedx: a true/false vector indicating which variables are fixed by
%    bound constraints
% 4. fixedx_value: the values of the variables fixed by bound constraints
% 5. nofreex: whether all variables are fixed by bound constraints
% 6. infeasible_bound: a true/false vector indicating which bound constraints
%    are infeasible
% 7. infeasible_lineq: a true/false vector indicating which linear inequality
%    constraints are infeasible (up to naive tests)
% 8. infeasible_leq: a true/false vector indicating which linear equality
%    constraints are infeasible (up to naive tests)
% 9. trivial_lineq
% 10. trivial_leq: a true/false vector indicating which linea equality
%     constraints are trivial (up to naive tests)
% 11. infeasible: whether the problem is infeasible (up to naive tests)
% 12. scaled: whether the problem is scaled
% 13. scaling_factor: vector of scaling factors
% 14. shift: vector of shifts
% 15. reduced: whether the problem is reduced (due to fixed variables)
% 16. raw_type: problem type before reduction
% 17. raw_dim: problem dimension before reduction
% 18. refined_type: problem type after reduction
% 19. refiend_dim: problem dimension after reduction
% 20. feasibility_problem: whether the problem is a feasibility problem
% 21. user_options_fields: the fields in the user-specified options
% 22. options: (refined) options for calling the solvers
% 23. warnings: warnings during the preprocessing/validation
probinfo = struct();

% Save the raw data (date before validation/preprocessing) in probinfo.
% The raw data can be useful when debugging. At the end of prepdfo, if
% we are not in debug mode, raw_data will be removed from probinfo.
% NOTE: Surely, here we are making copies of the data, which may take some
% time and space, matrices Aineq and Aeq being the major concern.
% However, fortunately, this package is not intended for large problems.
% It is designed for problems with at most ~1000 variables and several
% thousands of constraints, tens/hundreds of variables and tens/hundreds
% of constraints being typical. Therefore, making several copies (<10) of
% the data does not do much harm, especially when we solve problems with
% expensive (temporally or monetarily) function evaluations.
probinfo.raw_data = struct('objective', fun, 'x0', x0, 'Aineq', Aineq, 'bineq', bineq, ...
    'Aeq', Aeq, 'beq', beq, 'lb', lb, 'ub', ub, 'nonlcon', nonlcon, 'options', options);

% Validate and preprocess fun
[fun, probinfo.feasibility_problem, warnings] = pre_fun(invoker, fun, warnings);

% Validate and preprocess x0
[x0, warnings] = pre_x0(invoker, x0, warnings);
lenx0 = length(x0); % Within this file, for clarity, we denote length(x0) by lenx0 instead of n

% Validate and preprocess the bound constraints
% In addition, get the indices of infeasible bounds and 'fixed variables'
% such that ub-lb < 2eps (if any) and save the information in probinfo.
% If there is any infeasible bound, the problem is infeasible, and we define
% that there is no fixed variable.
[lb, ub, infeasible_bound, fixedx, fixedx_value, warnings] = pre_bcon(invoker, lb, ub, lenx0, warnings);
probinfo.infeasible_bound = infeasible_bound; % A vector of true/false
probinfo.fixedx = fixedx; % A vector of true/false
fixedx_value_save = fixedx_value; % Values of fixed variables
% Since fixedx_value may be revised in pre_lcon, we will record it in
% probinfo only after that. We save its current value in
% fixedx_value_save, which will be used when calculating the constraint
% violation at x0.

% Problem type before reduction
% This has to be done after preprocessing the bound constraints (because
% min(ub) and % max(lb) are evaluated) and before preprocessing the
% linear/nonlinear constraints (because these constraints will be
% reduced during the preprocessing). Note that Aineq, Aeq, and nonlcon will
% not be "evaluated" in problem_type, so there is no worry about the
% validity of them.
probinfo.raw_type = problem_type(Aineq, Aeq, lb, ub, nonlcon);

% Validate and preprocess the linear constraints
% 1. The constraints will be reduced if some but not all variables are
%    fixed by the bound constraints. See pre_lcon for why we do not
%    reduce the problem when all variables are fixed.
% 2. The 'trivial constraints' will be excluded (if any).
% 3. In addition, get the indices of infeasible and trivial constraints (if any)
%    and save the information in probinfo.
[Aineq, bineq, Aeq, beq, infeasible_lineq, trivial_lineq, infeasible_leq, trivial_leq, fixedx_value, warnings] = pre_lcon(invoker, x0, Aineq, bineq, Aeq, beq, lenx0, fixedx, fixedx_value, warnings);
probinfo.fixedx_value = fixedx_value; % Value of the fixed x entries; it is revised to the corresponding values of x0 in pre_lcon if infeasibility is detected
probinfo.infeasible_lineq = infeasible_lineq; % A vector of true/false
probinfo.trivial_lineq = trivial_lineq; % A vector of true/false
probinfo.infeasible_leq = infeasible_leq; % A vector of true/false
probinfo.trivial_leq = trivial_leq; % A vector of true/false

% Validate and preprocess the nonlinear constraints
% This should be done before evaluating probinfo.constrv_x0 or probinfo.constrv_fixedx.
% The constraints will be reduced if some but not all variables are fixed by the bound
% constraints. See pre_lcon for why we do not reduce the problem when all variables
% are fixed.
nonlcon = pre_nonlcon(invoker, nonlcon, fixedx, fixedx_value);

% Reduce fun, x0, lb, and ub if some but not all variables are fixed by
% the bound constraints. See pre_lcon for why we do not reduce the
% problem when all variables are fixed.
probinfo.raw_dim = lenx0; % Problem dimension before reduction
if any(fixedx) && any(~fixedx)
    freex = ~fixedx; % A vector of true/false indicating whether the variable is free or not
    fun = @(freex_value) fun(fullx(freex_value, fixedx_value, freex, fixedx)); % Objective funp after reduction
    x0 = x0(freex); % x0 after reduction
    lenx0 = length(x0);
    lb = lb(freex); % lb after reduction
    ub = ub(freex); % ub after reduction
end
probinfo.refined_type = problem_type(Aineq, Aeq, lb, ub, nonlcon); % Problem type after reduction
probinfo.refined_dim = length(x0); % Problem dimension after reduction
probinfo.reduced = any(fixedx) && any(~fixedx); % Whether the problem has been reduced

% After the preprocessing, the problem may turn out infeasible, or x may
% turn out fixed by the bounds
if ~any([probinfo.infeasible_lineq; probinfo.infeasible_leq; probinfo.infeasible_bound])
    probinfo.infeasible = false;
else % The problem turns out infeasible
    [probinfo.constrv_x0, probinfo.nlcineq_x0, probinfo.nlceq_x0] = constrv(x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon);
    % The constraint violation calculated by constrv does not include
    % the violation of x0 for the bounds corresponding to fixedx; the
    % corresponding values of x0 are in fixedx_value, while the values
    % of the bounds (lb and ub are the same up to eps) are in
    % fixedx_value_save. Thus the violation is abs(fixedx_value-fixedx_value_save).
    probinfo.constrv_x0 = max([probinfo.constrv_x0; abs(fixedx_value - fixedx_value_save)]);
    probinfo.infeasible = true;
end
if any(~fixedx)
    probinfo.nofreex = false;
else % x turns out fixed by the bound constraints
    [probinfo.constrv_fixedx, probinfo.nlcineq_fixedx, probinfo.nlceq_fixedx] = constrv(probinfo.fixedx_value, Aineq, bineq, Aeq, beq, lb, ub, nonlcon);
    probinfo.nofreex = true;
end

% Can the invoker handle the given problem?
% This should be done after the problem type has bee 'refined'.
if ~prob_solv_match(probinfo.refined_type, invoker)
    if strcmp(invoker, 'pdfon') || (nargin ~= 1)
        % Private/unexpected error
        error(sprintf('%s:InvalidProb', funname), ...
        '%s: UNEXPECTED ERROR: problem and solver do not match; this should not happen when %s is called by %s or the problem is not a structure.', funname, funname, invoker);
    else
        % Public/normal error
        error(sprintf('%s:InvalidProb', invoker), ...
        '%s: %s problem received; %s cannot solve it.', invoker, strrep(probinfo.refined_type, '-', ' '), invoker);
    end
end

% Validate and preprocess options, adopting default options if needed.
% This should be done after reducing the problem, because BOBYQA
% requires rhobeg <= min(ub-lb)/2.
% user_options_fields is a cell array that contains the names of all the
% user-defined options (even if the options turns out invalid). It will be
% needed if the user does not specify a solver or specifies a wrong solver.
% In such a scenario, we will select the solver later, and the options may
% have to be revised accordingly. We will raise a warning when revising
% an option that is in user_options_fields. No warning is needed if we
% are dealing with an option that is not in user_options_fields.
[options, probinfo.user_options_fields, warnings] = pre_options(invoker, options, lenx0, lb, ub, warnings);

% Revise x0 for bound and linearly constrained problems
% This is necessary for LINCOA, which accepts only feasible x0.
% Should we do this even if there are nonlinear constraints?
% For now, we do not, because doing so may dramatically increase the
% infeasibility of x0 with respect to the nonlinear constraints.
if ismember(probinfo.refined_type, {'bound-constrained', 'linearly-constrained'}) && ~probinfo.nofreex && ~probinfo.infeasible
    x0_old = x0;
    % Another possibility for bound-constrained problems:
    % xind = (x0 < lb) | (x0 > ub);
    % x0(xind) = (lb(xind) + ub(xind))/2;
    x0 = project(Aineq, bineq, Aeq, beq, lb, ub, x0);
    if norm(x0_old-x0) > eps*max(1, norm(x0_old)) && ~probinfo.feasibility_problem
        % No warning about revising x0 if the problem is a linear feasibility problem
        % Note that the linearity is guaranteed by THE OUTER IF.
        wid = sprintf('%s:ReviseX0', invoker);
        wmsg = sprintf('%s: x0 is revised to satisfy the constraints.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    end
end

% Scale the problem if necessary and if intended.
% x_before_scaling = scaling_factor.*x_after_scaling + shift
% This should be done after revising x0, which can affect the shift.
probinfo.scaled = false;
probinfo.scaling_factor = ones(size(x0));
probinfo.shift = zeros(size(x0));
if options.scale && ~probinfo.nofreex && ~probinfo.infeasible
    [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, scaling_factor, shift, ~, warnings] = scale_problem(invoker, fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, warnings);
    % Scale and shift the problem so that all the bounds become [-1, 1]
    % It is done only if all variables have both lower and upper bounds
    probinfo.scaled = true;
    probinfo.scaling_factor = scaling_factor;
    probinfo.shift = shift;
end

% Record the refined data (excluding options) after preprocessing
% This has to be done before select_solver, because probinfo.refined_data.lb
% and probinfo.refined_data.ub will be used for defining rhobeg if bobyqan is selected.
probinfo.refined_data = struct('objective', fun, 'x0', x0, 'Aineq', Aineq, 'bineq', bineq, ...
    'Aeq', Aeq, 'beq', beq, 'lb', lb, 'ub', ub, 'nonlcon', nonlcon);

% Select a solver if invoker='pdfon'; record the solver in options.solver.
% Some options will be revised accordingly, including npt, rhobeg, rhoend.
% Of course, if the user-defined options.solver is valid, we accept it.
if strcmp(invoker, 'pdfon')
    [options, warnings] = select_solver(invoker, options, probinfo, warnings);
end

if strcmpi(options.solver, 'bobyqan') && ~probinfo.nofreex && ~probinfo.infeasible && ~probinfo.feasibility_problem
% The Fortran code of BOBYQA will revise x0 so that the distance between
% x0 and the inactive bounds is at least rhobeg. We do it here in order
% to raise a warning when such a revision occurs. After this, the
% Fortran code will not revise x0 again. If the options.honour_x0 = true,
% then we keep x0 unchanged and revise rhobeg if necessary.
    [x0, options, warnings] = pre_rhobeg_x0(invoker, x0, lb, ub, probinfo.user_options_fields, options, warnings);
    probinfo.refined_data.x0 = x0;  % x0 may have been revised.
end

% Record the options in probinfo
% This has to be done after select_solver, because select_solver updates
% options.solver, and possibly options.npt and options.rhobeg.
% Also, pre_rhobeg_x0 may change options.rhobeg and options.rhoend.
probinfo.options = options;
% We do NOT record options in probinfo.refined_data, because we do not
% carry refined_data with us unless in debug mode or the problem is scaled.

if probinfo.feasibility_problem && ~strcmp(probinfo.refined_type, 'nonlinearly-constrained')
% When the problem is a linear feasibility problem, PDFO will return the
% current x0, which has been revised by project. The constraint violation
% at x0 is needed to set the output. Note that there is no nonlinear
% constraint in this case.
    probinfo.constrv_x0 = constrv(x0, Aineq, bineq, Aeq, beq, lb, ub, []);
end

probinfo.warnings = warnings; % Record the warnings in probinfo

if ~options.debug % Do not carry the raw data with us unless in debug mode.
    probinfo.raw_data = struct();
    % Set this field to empty instead of remove it, because postpdfo
    % requires this field to exist.
end

if ~options.debug && ~probinfo.scaled
    % The refined data is used only when the problem is scaled. It can
    % also be useful when debugging.
    probinfo.refined_data = struct();
    % Set this field to empty instead of remove it, because postpdfo
    % requires this field to exist.
end

% prepdfo ends
return

%%%%%%%%%%%%%%%%%%%%%%%% Function for problem decoding %%%%%%%%%%%%%%%%%
function [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options, warnings] = decode_problem(invoker, problem, warnings)
% Read the fields of the 'problem' structure but do not validate them.
% The decoded problem will be sent to the prepdfo function for validation.
% NOTE: We treat field names case-sensitively.

% Possible invokers
invoker_list = {'pdfon', 'uobyqan', 'newuoan', 'bobyqan', 'lincoan', 'cobylan'};

callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if ~ismember(invoker, invoker_list)
    % invoker affects the behavior of this function, so we check invoker
    % again, even though it should have been checked in function prepdfo
    % Private/unexpcted error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s serves only %s.', funname, funname, mystrjoin(invoker_list, ', '));
end

if ~isa(problem, 'struct')
    % Public/normal error
    error(sprintf('%s:InvalidProb', invoker), '%s: the unique input is not a problem-defining structure.', invoker);
end

% Which fields are specified?
problem = rmempty(problem); % Remove empty fields
problem_fields = fieldnames(problem);

% Are the obligatory field(s) present?
obligatory_fields = {'x0'}; % There is only 1 obligatory field
missing_fields = setdiff(obligatory_fields, problem_fields);
if ~isempty(missing_fields)
    % Public/normal error
    error(sprintf('%s:InvalidProb', invoker), ...
    '%s: PROBLEM misses the %s field(s).', invoker, mystrjoin(missing_fields, ', '));
end
x0 = problem.x0;

if isfield(problem, 'objective')
    fun = problem.objective;
else % There is no objective; this is a feasibility problem
    fun = []; % We use [] to signify that fun is not specified. pre_fun will replace [] by by @(x)0
end

% Are there unknown fields?
known_fields = {'objective', 'x0', 'Aineq', 'bineq', 'Aeq', 'beq', 'lb', 'ub', 'nonlcon', 'options', 'solver'};
% 1. When invoker is in {uobyqan, ..., cobylan}, we will not complain that
%    a solver is specified unless invoker~=solver. See function pre_options.
% 2. When invoker is in {uobyqan, ..., cobylan}, if the problem turns out
%    unsolvable for the invoker, then we will raise an error in prepdfo.
%    We do not do it here because the problem has not been validated/preprocessed
%    yet. Maybe some constraints are trivial and hence can be removed
%    (e.g., bineq=inf, lb=-inf), which can change the problem type.

unknown_fields = setdiff(problem_fields, known_fields);
problem = rmfield(problem, unknown_fields);  % Remove the unknown fields

if ~isempty(unknown_fields)
    wid = sprintf('%s:UnknownProbField', invoker);
    if length(unknown_fields) == 1
        wmsg = sprintf('%s: problem with an unknown field %s; it is ignored.', invoker, mystrjoin(unknown_fields, ', '));
    else
        wmsg = sprintf('%s: problem with unknown fields %s; they are ignored.', invoker, mystrjoin(unknown_fields, ', '));
    end
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Read the fields of problem. They will be validated in function predfo
Aineq = [];
bineq = [];
Aeq = [];
beq = [];
lb= [];
ub = [];
nonlcon = [];
options = struct();
if isfield(problem,'Aineq')
    Aineq = problem.Aineq;
end
if isfield(problem,'bineq')
    bineq = problem.bineq;
end
if isfield(problem,'Aeq')
    Aeq = problem.Aeq;
end
if isfield(problem,'beq')
    beq = problem.beq;
end
if isfield(problem,'lb')
    lb = problem.lb;
end
if isfield(problem,'ub')
    ub = problem.ub;
end
if isfield(problem,'nonlcon')
    nonlcon = problem.nonlcon;
end
if isfield(problem,'options')
    options = problem.options;
end
if isfield(problem,'solver')
    options.solver = problem.solver;
    % After last step, options.solver = problem.options.solver;
    % after this step, if problem.solver is defined and nonempty,
    % then options.solver = problem.solver.
end
return

%%%%%%%%%%%%%%%%%%%%%%%% Function for fun preprocessing %%%%%%%%%%%%%%%%%
function [fun, feasibility_problem, warnings] = pre_fun(invoker, fun, warnings)
if ~(isempty(fun) || isa(fun, 'char') || isa(fun, 'string') || isa(fun, 'function_handle'))
    % Public/normal error
    error(sprintf('%s:InvalidFun', invoker), ...
        '%s: FUN should be a function handle or a function name.', invoker);
end
feasibility_problem = false; % Is this a feasibility problem?
if isempty(fun)
    fun = @(x) 0; % No objective function
    feasibility_problem = true; % This is a feasibility problem
    wid = sprintf('%s:NoObjective', invoker);
    wmsg = sprintf('%s: there is no objective function. A feasibility problem will be solved.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
elseif isa(fun, 'char') || isa(fun, 'string')
    fun = str2func(fun);
    % Work with function handles instread of function names to avoid using 'feval'
end
if ~exist('OCTAVE_VERSION', 'builtin')
    % Check whether fun has at least 1 output.
    % nargout(fun) = #outputs in the definition of fun.
    % If fun includes varargout in definition, nargout(fun) = -#outputs.
    % Octave does not support nargout for built-in function (as of 2019-08-16)!
    try
    % If fun is not a properly defined function, then nargout
    % can encounter an error. Wrap the error as a public error.
        nout = nargout(fun);
    catch exception
        % Public/normal error
        % Note that the identifier of a public error should start with 'invoker:'
        error(sprintf('%s:InvalidFun', invoker), '%s: %s', invoker, exception.message);
    end
    if (nout == 0)
        % Public/normal error
        error(sprintf('%s:InvalidFun', invoker), ...
        '%s: FUN has no output; it should return the objective function value.', invoker);
    end
end
fun = @(x) evalobj(invoker, fun, x);
return

function f = evalobj(invoker, fun, x)
f = fun(x);
if ~isnumeric(f) || numel(f) ~= 1
    % Public/normal error
    error(sprintf('%s:ObjectiveNotScalar', invoker), '%s: objective function should return a scalar value.', invoker);
end
f = double(real(f)); % Some functions like 'asin' can return complex values even when it is not intended
% Use extreme barrier to cope with 'hidden constraints'
hugefun = gethuge('fun');
if (f ~= f) || (f > hugefun)
    f = hugefun;
end
return

%%%%%%%%%%%%%%%%%%%%%%%% Function for x0 preprocessing %%%%%%%%%%%%%%%%%
function [x0, warnings] = pre_x0(invoker, x0, warnings)
[isrv, lenx0]  = isrealvector(x0);
if ~(isrv && (lenx0 > 0))
    % Public/normal error
    error(sprintf('%s:InvalidX0', invoker), '%s: X0 should be a real vector/scalar.', invoker);
end
x0 = double(x0(:));
abnormal_x0 = isnan(x0) | (abs(x0) >= inf);
if any(abnormal_x0)
    x0(abnormal_x0) = 0;
    wid = sprintf('%s:AbnormalX0', invoker);
    wmsg = sprintf('%s: X0 contains NaN or inifinite values; they are replaced by 0.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end
return

%%%%%%%%%%%%%%%%% Function for bound constraint preprocessing %%%%%%%%%%
function [lb, ub, infeasible_bound, fixedx, fixedx_value, warnings] = pre_bcon(invoker, lb, ub, lenx0, warnings)
% Lower bounds (lb)
[isrvlb, lenlb] = isrealvector(lb);
if ~(isrvlb && (lenlb == lenx0 || lenlb == 0))
    % Public/normal error
    error(sprintf('%s:InvalidBound', invoker), ...
    '%s: lb should be a real vector and length(lb)=length(x0) unless lb=[].', invoker);
end
if (lenlb == 0)
    lb = -inf(lenx0,1); % After pre_bcon, length(lb) = length(x0)
end
lb = double(lb(:));
if any(isnan(lb))
    lb(isnan(lb)) = -inf; % Replace the NaN in lb by -inf
    wid = sprintf('%s:NaNInLB', invoker);
    wmsg = sprintf('%s: LB contains NaN; it is replaced by -inf.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Upper bounds (ub)
[isrvub, lenub] = isrealvector(ub);
if ~(isrvub && (lenub == lenx0 || lenub == 0))
    % Public/normal error
    error(sprintf('%s:InvalidBound', invoker), ...
    '%s: ub should be a real vector and length(ub)=length(x0) unless ub=[].', invoker);
end
if (lenub == 0)
    ub = inf(lenx0,1); % After pre_bcon, length(ub) = length(x0)
end
ub = double(ub(:));
if any(isnan(ub))
    ub(isnan(ub)) = inf; % Replace the NaN in ub by inf
    wid = sprintf('%s:NaNInUB', invoker);
    wmsg = sprintf('%s: UB contains NaN; it is replaced by inf.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

infeasible_bound = (lb > ub) | (lb == inf) | (ub == -inf); % A vector of true/false
if any(infeasible_bound)
    fixedx = false(lenx0, 1);
    fixedx_value = [];
else
    fixedx = (abs(lb - ub) < 2*eps);
    fixedx_value = (lb(fixedx)+ub(fixedx))/2;
end
return

%%%%%%%%%%%%%%%%% Function for linear constraint preprocessing %%%%%%%%%%
function [Aineq, bineq, Aeq, beq, infeasible_lineq, trivial_lineq, infeasible_leq, trivial_leq, fixedx_value, warnings] = pre_lcon(invoker, x0, Aineq, bineq, Aeq, beq, lenx0, fixedx, fixedx_value, warnings)

freex = ~fixedx; % A vector of true/false indicating whether the variable is free or not

% inequalities: Aineq*x <= bineq
[isrm, mA, nA] = isrealmatrix(Aineq);
[isrc, lenb] = isrealcolumn(bineq);
if ~(isrm && isrc && (mA == lenb) && (nA == lenx0 || nA == 0))
    % Public/normal error
    error(sprintf('%s:InvalidLinIneq', invoker), ...
    '%s: Aineq should be a real matrix, bineq should be a real column, and size(Aineq)=[length(bineq), length(X0)] unless Aineq=bineq=[].', invoker);
end
if any(isnan(bineq))
    bineq(isnan(bineq)) = inf; % Replace the NaN in bineq by inf
    wid = sprintf('%s:NaNInbineq', invoker);
    wmsg = sprintf('%s: bineq contains NaN; it is replaced by inf.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end
lineq_reduced = false; % Whether linear inequality constraints are reduced
if ~isempty(Aineq) && any(fixedx) && any(~fixedx)
    % Reduce the linear inequality constraints if some but not all variables
    % are fixed by the bound constraints. This has to be done before
    % detecting the "zero constraints" (i.e., constraints with zero
    % gradients), because nonzero constraints may become zero after reduction.
    Aineq_fixed = Aineq(:, fixedx); % Aineq_fixed and bineq_save will be used when revising fixedx_value
    bineq_save = bineq;
    bineq = bineq - Aineq_fixed * fixedx_value;
    Aineq = Aineq(:, freex);
    lineq_reduced = true;
    % Note that we should NOT reduced the problem if all variables are
    % fixed. Otherwise, Aineq would be [], and then bineq will be set to
    % [] in the end. In this way, we lose completely the information in
    % linear constraints. Consequently, we cannot evaluate the constraint
    % violation correctly when needed.
end
if isempty(Aineq)
    infeasible_lineq = [];
    trivial_lineq = [];
else
    Aineq = double(Aineq);
    bineq = double(bineq);
    rownorm1 = sum(abs(Aineq), 2);
    zero_ineq = (rownorm1 == 0);
    infeasible_zero_ineq = (rownorm1 == 0) & (bineq < 0);
    trivial_zero_ineq = (rownorm1 == 0) & (bineq >= 0);
    rownorm1(zero_ineq) = 1;
    infeasible_lineq = (bineq./rownorm1 == -inf) | infeasible_zero_ineq | isnan(rownorm1); % A vector of true/false
    trivial_lineq = (bineq./rownorm1 == inf) | trivial_zero_ineq;
    Aineq = Aineq(~trivial_lineq, :); % Remove the trivial linear inequalities
    bineq = bineq(~trivial_lineq);
end

% equalities: Aeq*x == beq
[isrm, mA, nA] = isrealmatrix(Aeq);
[isrc, lenb] = isrealcolumn(beq);
if ~(isrm && isrc && (mA == lenb) && (nA == lenx0 || nA == 0))
    % Public/normal error
    error(sprintf('%s:InvalidLinEq', invoker), ...
    '%s: Aeq should be a real matrix, beq should be a real column, and size(Aeq)=[length(beq), length(X0)] unless Aeq=beq=[].', invoker);
end
% Are there equality constraints whose both sides contain NaN?
% This should be detected before reducing the constraints;
% when reducing the constraints, the NaN on the left-hand side will lead
% to NaN on the right-hand side.
if isempty(Aeq)
    nan_eq = [];
else
    nan_eq = isnan(sum(abs(Aeq), 2)) & isnan(beq); % In MATLAB 2014a, this may lead to inconsistent sizes when Aeq is empty; in MATLAB 2018a, it is fine
end
if any(nan_eq)
    wid = sprintf('%s:NaNEquality', invoker);
    wmsg = sprintf('%s: there are equality constraints whose both sides contain NaN; such constraints are removed.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end
leq_reduced = false; % Whether linear equality constraints are reduced
if ~isempty(Aeq) && any(fixedx) && any(~fixedx)
    % Reduce the linear equality constraints if some but not all variables
    % are fixed by the bound constraints. This has to be done before
    % detecting the "zero constraints" (i.e., constraints with zero
    % gradients), because nonzero constraints may become zero after reduction.
    Aeq_fixed = Aeq(:, fixedx); % Aeq_fixed and beq_save may be used when revising fixedx_value
    beq_save = beq;
    beq = beq - Aeq_fixed * fixedx_value;
    Aeq = Aeq(:, freex);
    leq_reduced = true;
    % Note that we should NOT reduced the problem if all variables are
    % fixed. Otherwise, Aeq would be [], and then beq will be set to
    % [] in the end. In this way, we lose completely the information in
    % linear constraints. Consequently, we cannot evaluate the constraint
    % violation correctly when needed.
end
if isempty(Aeq)
    infeasible_leq = [];
    trivial_leq = [];
else
    Aeq = double(Aeq);
    beq = double(beq);
    rownorm1 = sum(abs(Aeq), 2);
    zero_eq = (rownorm1 == 0);
    infeasible_zero_eq = (rownorm1 == 0) & (beq ~= 0);
    trivial_zero_eq = (rownorm1 == 0) & (beq == 0);
    rownorm1(zero_eq) = 1;
    infeasible_leq = (abs(beq./rownorm1) == inf) | infeasible_zero_eq | ((isnan(rownorm1) | isnan(beq)) & ~nan_eq); % A vector of true/false
    trivial_leq = trivial_zero_eq | nan_eq;
    Aeq = Aeq(~trivial_leq, :); % Remove trivial linear equalities
    beq = beq(~trivial_leq);
end

% In case of infeasibility, revise fixedx_value so that x0 will be returned
if (any(infeasible_lineq) || any(infeasible_leq)) && any(fixedx) && any(~fixedx)
    fixedx_value = x0(fixedx);
    % We have to revise bineq and beq so that the constraint violation can
    % be correctly calculated.
    if lineq_reduced
        bineq = bineq_save - Aineq_fixed * fixedx_value;
        bineq = bineq(~trivial_lineq);
    end
    if leq_reduced
        beq = beq_save - Aeq_fixed * fixedx_value;
        beq = beq(~trivial_leq);
    end
end

% We uniformly use [] to represent empty numerical matrices/vectors;
% its size is 0x0. Changing this may cause matrix dimension inconsistency.
if isempty(Aeq)
    Aeq = [];
    beq = [];
end
if isempty(Aineq)
    Aineq = [];
    bineq = [];
end
return

%%%%%%%%%%%%%%%%% Function for nonlinear constraint preprocessing %%%%%%%%%%
function nonlcon = pre_nonlcon(invoker, nonlcon, fixedx, fixedx_value)
if ~(isempty(nonlcon) || isa(nonlcon, 'function_handle') || isa(nonlcon, 'char') || isa(nonlcon, 'string'))
    % Public/normal error
    error(sprintf('%s:InvalidCon', invoker), ...
    '%s: nonlcon should be a function handle or a function name.', invoker);
end
if isempty(nonlcon)
    nonlcon = []; % We use [] to signify that nonlcon is not specified; its size is 0x0
else
    if isa(nonlcon, 'char') || isa(nonlcon, 'string')
        nonlcon = str2func(nonlcon);
        % work with function handles instead of function names to avoid using 'feval'
    end
    if ~exist('OCTAVE_VERSION', 'builtin')
        % Check whether nonlcon has at least 2 outputs.
        % nargout(fun) = #outputs in the definition of fun.
        % If fun includes varargout in definition, nargout(fun) = -#outputs.
        % Octave does not support nargout for built-in function (as of 2019-08-16)!
        try
        % If nonlcon is not a properly defined function, then nargout
        % can encounter an error. Wrap the error as a public error.
            nout = nargout(nonlcon);
        catch exception
            % Public/normal error
            % Note that the identifier of a public error should start with 'invoker:'
            error(sprintf('%s:InvalidCon', invoker), '%s: %s', invoker, exception.message);
        end
        if (nout == 0) || (nout == 1)
            % Public/normal error
            error(sprintf('%s:InvalidCon', invoker), ...
            '%s: nonlcon has too few outputs; it should return [cineq, ceq], the constraints being cineq(x)<=0, ceq(x)=0.', invoker);
        end
    end
    if any(fixedx) && any(~fixedx)
        % Reduce the nonlinear constraints if some but not all variables are
        % fixed by the bound constraints. Note that we do not reduce the
        % problem when all variables are fixed. See pre_lcon for the reason.
        freex = ~fixedx; % A vector of true/false indicating whether the variable is free or not
        nonlcon = @(freex_value) nonlcon(fullx(freex_value, fixedx_value, freex, fixedx));
    end
    nonlcon = @(x) evalcon(invoker, nonlcon, x);
end
return

function [cineq, ceq] = evalcon(invoker, nonlcon, x)
[cineq, ceq] = nonlcon(x);
if ~(isempty(cineq) || isnumeric(cineq)) || ~(isempty(ceq) || isnumeric(ceq))
    % Public/normal error
    error(sprintf('%s:ConstrNotNumeric', invoker), '%s: constraint function should return two numeric vectors.', invoker);
end
cineq = double(real(cineq(:))); % Some functions like 'asin' can return complex values even when it is not intended
ceq = double(real(ceq(:)));
% Use extreme barrier to cope with 'hidden constraints'
hugecon = gethuge('con');
cineq(cineq ~= cineq) = hugecon;
cineq(cineq > hugecon) = hugecon;
ceq(ceq ~= ceq) = hugecon;
ceq(ceq > hugecon) = hugecon;
ceq(ceq < -hugecon) = -hugecon;

% This part is NOT extreme barrier. We replace extremely negative values of
% cineq (which leads to no constraint violation) by -hugecon. Otherwise,
% NaN or Inf may occur in the interpolation models.
cineq(cineq < -hugecon) = -hugecon;
return

%%%%%%%%%%%%%%%%% Function fullx used when reducing the problem %%%%%%%%
function x = fullx(freex_value, fixedx_value, freex, fixedx)
x = NaN(length(freex_value)+length(fixedx_value), 1);
x(freex) = freex_value;
x(fixedx) = fixedx_value;
return

%%%%%%%%%%%%%%%%% Function for option preprocessing %%%%%%%%%%
function [options, user_options_fields, warnings] = pre_options(invoker, options, lenx0, lb, ub, warnings)

% NOTE: We treat field names case-sensitively.

% Possible solvers
solver_list = {'uobyqan', 'newuoan', 'bobyqan', 'lincoan', 'cobylan'};
% We may add other solvers in the future!
% If a new solver is included, we should do the following.
% 0. Include it into the invoker_list (in this and other functions).
% 1. What options does it expect? Set known_fields accordingly.
% 2. Set default options accordingly.
% 3. Check other functions (especially decode_problem, whose behavior
%    depends on the invoker/solver. See known_fields there).

% Possible invokers
invoker_list = ['pdfon', solver_list];

callstack = dbstack;
funname = callstack(1).name;
% invoker affects the behavior of this function, so we check invoker
% again, even though it should have been checked in function prepdfo
if ~ismember(invoker, invoker_list)
    % Private/unexpcted error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s serves only %s.', funname, funname, mystrjoin(invoker_list, ', '));
end

% Default values of the options.
% npt = ! LATER ! % The default npt depends on solver and will be set later in this function
maxfun = 500*lenx0;
rhobeg = 1; % The default rhobeg and rhoend will be revised if solver = 'bobyqan'
rhoend = 1e-6;
ftarget = -inf;
classical = false; % Call the classical Powell code? Classical mode recommended only for research purpose
fortran = true; % Call the Fortran code?
scale = false; % Scale the problem according to bounds? Scale only if the bounds reflect well the scale of the problem
scale = (scale && max(ub-lb)<inf); % ! NEVER remove this ! Scale only if all variables are with finite lower and upper bounds
honour_x0 = false; % Respect the user-defined x0? Needed by BOBYQA
iprint = 0;
quiet = true;
debugflag = false; % Do not use 'debug' as the name, which is a MATLAB function
chkfunval = false;
ctol = eps; % Tolerance for constraint violation; a point with a constraint violation at most ctol is considered feasible
output_xhist = false; % Output the history of x?
output_nlchist = false; % Output the history of the nonlinear constraints?
maxfilt = 2000; % Length of the filter used for selecting the returned x in constrained problems
min_maxfilt = 200; % The smallest value of maxfilt; if maxfilt is too small, the returned x may not be the best one visited

if ~(isa(options, 'struct') || isempty(options))
    % Public/normal error
    error(sprintf('%s:InvalidOptions', invoker), '%s: OPTIONS should be a structure.', invoker);
end

% Which fields are specified?
options = rmempty(options); % Remove empty fields
options_fields = fieldnames(options);
% The list of fields in options  will be returned and used elsewhere. We
% save it right now in case we "intellegently" change options_fields
% after this line in future versions.
user_options_fields = options_fields;

% Validate options.solver
% We need to know what is the solver in order to decide which fields
% are 'known' (e.g., expected), and also to set npt, rhobeg, rhoend.
% We do the following:
% 1. If invoker='pdfon':
% 1.1 If no solver is specified or solver='pdfon', we do not complain
% and set options.solver=solver='', i.e., an empty char array;
% 1.2 Else if solver is not in solver_list, we warn about 'unknown solver'
% and set options.solver=solver='', i.e., an empty char array;
% 1.3 Else, we set solver=options.solver.
% 2. If invoker is in solver_list:
% 2.1 If options.solver exists but options.solver~=invoker, we warn
% about 'invalid solver' and set options.solver=solver=invoker;
% 2.2 Else, we do not complain and set options.solver=solver=invoker.
% In this way, options.solver and solver either end up with a member of
% solver_list or ''. The second case is possible only if invoker='pdfon',
% and solver will be selected later.
if isfield(options, 'solver') && ~isa(options.solver, 'char') && ~isa(options.solver, 'string')
    options.solver = 'UNKNOWN_SOLVER';
    % We have to change options.solver to a char/string so that we can use strcmpi
    % We do not need to worry about the case where solver is empty, because
    % all the empty fields have been removed from options.
end
if strcmp(invoker, 'pdfon')
    % We se the default value of solver to '', an empty char array.
    % 1. DO NOT change this default value! It will affect known_fields
    % and select_solver.
    % 2. DO NOT use [], which is an empty double array and may cause some
    % functions (e.g., ismember) to complain about incompatible types.
    solver = '';
    if isfield(options, 'solver')
        if any(strcmpi(options.solver, solver_list))
            solver = lower(options.solver);
        elseif ~strcmpi(options.solver, 'pdfon')
        % We should not complain about 'unknown solver' if invoker=options.solver='pdfon'
            wid = sprintf('%s:UnknownSolver', invoker);
            wmsg = sprintf('%s: unknown solver specified; %s will select one automatically.', invoker, invoker);
            warning(wid, '%s', wmsg);
            warnings = [warnings, wmsg];
        end
    end
else % invoker is in {'uobyqan', ..., 'cobylan'}
    if isfield(options, 'solver') && ~strcmpi(options.solver, invoker)
        wid = sprintf('%s:InvalidSolver', invoker);
        wmsg = sprintf('%s: a solver different from %s is specified; it is ignored.', invoker, invoker);
        % Do not display the value of solver in last message, because it
        % can be 'unknow_solver'.
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    end
    solver = invoker;
end
options.solver = solver; % Record solver in options.solver; will be used in postpdfo
% When the invoker is pdfo, options.solver=solver='' unless the user defines
% an options.solver in solver_list. Here, '' is an empty char array to signify
% that the solver is yet to decide.

% Check unknown fields according to solver
% solver is '' if it has not been decided yet; in that case, we suppose (for
% simplicity) that all possible fields are known.
known_fields = {'iprint', 'maxfun', 'rhobeg', 'rhoend', 'ftarget', 'classical', 'quiet', 'debug', 'chkfunval', 'solver', 'maxhist', 'output_xhist', 'fortran'};
if ~isfield(options, 'classical') || (islogicalscalar(options.classical) && ~options.classical)
    known_fields = [known_fields, 'eta1', 'eta2', 'gamma1', 'gamma2'];
end
if isempty(solver) || any(strcmpi(solver, {'newuoan', 'bobyqan', 'lincoan'}))
    known_fields = [known_fields, 'npt'];
end
if isempty(solver) || any(strcmpi(solver, {'bobyqan', 'lincoan', 'cobylan'}))
    known_fields = [known_fields, 'scale'];
end
if isempty(solver) || strcmpi(solver, 'bobyqan')
    known_fields = [known_fields, 'honour_x0'];
end
if isempty(solver) || any(strcmpi(solver, {'lincoan', 'cobylan'}))
    known_fields = [known_fields, {'ctol', 'maxfilt'}];
end
if isempty(solver) || strcmpi(solver, 'cobylan')
    known_fields = [known_fields, 'output_nlchist'];
end
unknown_fields = setdiff(options_fields, known_fields);
options = rmfield(options, unknown_fields);  % Remove the unknown fields
% If we do not removed unknown fields, we may still complain later if an
% unknown field is not properly set (e.g., options.npt is not a number)
% even though we have declared that this field will be ignored.
if ~isempty(unknown_fields)
    wid = sprintf('%s:UnknownOption', invoker);
    if length(unknown_fields) == 1
        wmsg = sprintf('%s: unknown option %s; it is ignored.', invoker, mystrjoin(unknown_fields, ', '));
    else
        wmsg = sprintf('%s: unknown options %s; they are ignored.', invoker, mystrjoin(unknown_fields, ', '));
    end
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Set default npt according to solver
% If solver='' (empty char array), then invoker must be pdfo, and a solver
% will be selected later; when the solver is chosen, a valid npt will be defined.
% Note we have to take maxfun into consideration when selecting the solver,
% because npt < maxfun-1 is needed! See function select_solver for details.
if isempty(solver)
    npt = NaN; % The real npt will be (and should be) set when solver is selected
else
    switch lower(solver)
    case {'newuoan', 'bobyqan', 'lincoan'}
        npt = 2*lenx0 + 1;
    case {'uobyqan'}
        npt = (lenx0+1)*(lenx0+2)/2;
    case {'cobylan'}
        npt = lenx0+1;
		% uobyqan and cobylan do not need npt an option, but we need npt to validate/set maxfun
    end
end

% Validate options.scale
% We need the value of options.scale to revise the default rhobeg
validated = false;
if isfield(options, 'scale')
    if ~islogicalscalar(options.scale)
        wid = sprintf('%s:InvalidScaleFlag', invoker);
        wmsg = sprintf('%s: invalid scale flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(scale));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif options.scale && max(ub-lb) >= inf
        wid = sprintf('%s:ProblemCannotBeScaled', invoker);
        wmsg = sprintf('%s: problem cannot be scaled because not all variables have both lower and upper bounds.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
        options.scale = false; % options.scale must be set to false in this case
        validated = true;
    else
        validated = true;
    end
end
if ~validated % options.scale has not got a valid value yet
    options.scale = scale;
end
options.scale = logical(options.scale);

% Revise default rhobeg and rhoend according to options.scale and solver
if options.scale
    rhobeg = 0.5;  % This value cannot be bigger than 1. Otherwise, BOBYQA will complain.
    rhoend = 1e-6;
end
if strcmpi(solver, 'bobyqan') && ~options.scale
    rhobeg = max(eps, min(rhobeg, min(ub-lb)/4));
    rhoend = max(eps, min(0.1*rhobeg, rhoend));
end


% Validate the user-specified options; adopt the default values if needed

% Validate options.npt
% There are the following possibilities.
% 1. The user specifies options.npt
% 1.1. The solver is yet to decide (solver=''): we keep options.npt if it is
% a positive integer; otherwise, raise a warning and set options.npt to NaN;
% 1.2. The user has chosen a valid solver: we keep options.npt if it is
% compatible with the solver; otherwise, raise a warning and set options.npt
% to the default value according to the solver.
% 2. The user does not specify options.npt
% 1.1. The solver is yet to decide (solver=''): we set options.npt to NaN.
% 1.2. The user has chosen a valid solver: we set options.npt to the default
% value accoring to the solver.
% After this process, options.npt is either a positive integer (compatible
% with options.solver if it is specified by the user) or NaN (only if the
% user does not specify a valid solver while options.npt is either unspecified
% or not a positive integer).
validated = false;
if isfield(options, 'npt')
    if isempty(solver) && (~isintegerscalar(options.npt) || options.npt < 1 || isnan(options.npt))
        wid = sprintf('%s:InvalidNpt', invoker);
        wmsg = sprintf('%s: invalid npt. It should be a positive integer.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif any(strcmpi(solver, {'newuoan', 'bobyqan', 'lincoan'})) && (~isintegerscalar(options.npt) || isnan(options.npt) || options.npt < lenx0+2 || options.npt > (lenx0+1)*(lenx0+2)/2)
        % newuoan, bobyqan and lincoan requires n+2<=npt<=(n+1)*)(n+2)/2;
        % uobyqan and cobylan do not use npt.
        wid = sprintf('%s:InvalidNpt', invoker);
        wmsg = sprintf('%s: invalid npt; %s requires it to be an integer and n+2 <= npt <= (n+1)*(n+2)/2; it is set to 2n+1.', invoker, solver);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated  % options.npt has not got a valid value yet
    options.npt = npt;
    % When solver='' (empty char array), the default npt is NaN.
    % For uobyqan and cobylan, we also adopt the 'default npt' defined above,
    % although it will NOT be used by the solver
end
options.npt = double(options.npt);
% Although npt and maxfun are integers logically, they have to be
% passed to the mexified code as double variables. In mex, data is
% passed by pointers, but there are only very limited functions that
% can read an integer value from a pointer or write an interger
% value to a pointer (mxCopyPtrToInteger1, mxCopyInteger1ToPtr,
% mxCopyPtrToInteger2, mxCopyInteger2ToPtr, mxCopyPtrToInteger4,
% mxCopyInteger4ToPtr; no function for integer*8). This makes it
% impossible to pass integer data properly unless we know the kind
% of the integer. Therefore, in general, it is recommended to pass
% integers as double variables and then cast them back to integers
% when needed.
% Indeed, in matlab, even if we define npt = 1000,
% the class of npt is double! To get an integer npt, we would
% have to define npt = int32(1000) or npt = int64(1000)!

% Validate options.maxfun
validated = false;
if isfield(options, 'maxfun')
    if ~isintegerscalar(options.maxfun) || options.maxfun <= 0 || isnan(options.maxfun) || options.maxfun == inf
        % Here, we do not revise excessively large maxfun (e.g., maxfun = 10^100),
        % which should be handled by each solver case by case.
        wid = sprintf('%s:InvalidMaxfun', invoker);
        wmsg = sprintf('%s: invalid maxfun; it should be a positive integer; it is set to %d.', invoker, maxfun);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif isempty(solver) && options.maxfun <= lenx0+1  % Here, options.maxfun cannot be NaN. No worry about the comparison.
        options.maxfun = lenx0+2; % Here we take lenx0+2 (the smallest possible value for npt)
        validated = true; %!!! % Set validated=true so that options.maxfun will not be set to the default value later
        wid = sprintf('%s:InvalidMaxfun', invoker);
        wmsg = sprintf('%s: invalid maxfun; it should be a positive integer at least n+2; it is set to n+2.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif ~isempty(solver) && options.maxfun <= options.npt  % Here, options.maxfun or options.npt cannot be NaN. No worry about the comparison.
        options.maxfun = options.npt+1; % Here we take npt+1 instead of the default maxfun
        validated = true; %!!! % Set validated=true so that options.maxfun will not be set to the default value later
        wid =  sprintf('%s:InvalidMaxfun', invoker);
        switch lower(solver) % The warning message depends on solver
        case {'newuoan', 'lincoan', 'bobyqan'}
            wmsg = sprintf('%s: invalid maxfun; %s requires maxfun > npt; it is set to npt+1.', invoker, solver);
        case 'uobyqan'
            wmsg = sprintf('%s: invalid maxfun; %s requires maxfun > (n+1)*(n+2)/2; it is set to (n+1)*(n+2)/2+1.', invoker, solver);
        case 'cobylan'
            wmsg = sprintf('%s: invalid maxfun; %s requires maxfun > n+1; it is set to n+2.', invoker, solver);
        end
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.maxfun has not got a valid value yet
    options.maxfun = max(maxfun, options.npt+1);
end
options.maxfun = double(options.maxfun); % maxfun will be passed as a double
% One can check that options.maxfun >= n+2;

% Validate options.rhobeg
% NOTE: if the problem is to be scaled, then options.rhobeg and options.rhoend
% will be used as the intial and final trust-region radii for the scaled problem.
validated = false;
if isfield(options, 'rhobeg')
    if ~isrealscalar(options.rhobeg) || options.rhobeg <= 0 || isnan(options.rhobeg) || options.rhobeg == inf
        wid = sprintf('%s:InvalidRhobeg', invoker);
        wmsg = sprintf('%s: invalid rhobeg; it should be a positive number; it is set to max(%f, rhoend).', invoker, rhobeg);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif strcmpi(solver, 'bobyqan')  % Validate options.rhobeg for bobyqan
        if options.scale && options.rhobeg > 1  % This case cannot be combined with the next case, as ub and lb are NOT scaled yet in prepdfo
            wid = sprintf('%s:InvalidRhobeg', invoker);
            wmsg = sprintf('%s: invalid rhobeg; %s requires rhobeg <= 1 when the problem is scaled; it is set to 0.5.', invoker, solver);
            warning(wid, '%s', wmsg);
            warnings = [warnings, wmsg];
            options.rhobeg = 0.5;
        elseif ~options.scale && options.rhobeg > min(ub-lb)/2
            wid = sprintf('%s:InvalidRhobeg', invoker);
            wmsg = sprintf('%s: invalid rhobeg; %s requires rhobeg <= min(ub-lb)/2; it is set to min(ub-lb)/4.', invoker, solver);
            warning(wid, '%s', wmsg);
            warnings = [warnings, wmsg];
            options.rhobeg = min(ub-lb)/4; % Here we do not take the default rhobeg
        end
        validated = true; %!!! % Set validated=true so that options.rhobeg will not be set to the default value later
    else
        validated = true;
    end
end
if ~validated % options.rhobeg has not got a valid value yet
    if isfield(options, 'rhoend') && isrealscalar(options.rhoend) && options.rhoend >=0 && options.rhoend < inf
        options.rhobeg = max(rhobeg, 10*options.rhoend);
    else
        options.rhobeg = rhobeg;
    end
end
options.rhobeg = double(max(options.rhobeg, eps));

% Validate options.rhoend
validated = false;
if isfield(options, 'rhoend')
    if ~isrealscalar(options.rhoend) || options.rhoend > options.rhobeg || isnan(options.rhoend)
        wid = sprintf('%s:InvalidRhoend', invoker);
        wmsg = sprintf('%s: invalid rhoend; we should have rhobeg >= rhoend > 0; it is set to min(0.1*rhobeg, %f).', invoker, rhoend);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.rhoend has not got a valid value yet
    options.rhoend = min(0.1*options.rhobeg, rhoend);
end
options.rhoend = double(max(options.rhoend, eps));

% Validate options.ftarget
validated = false;
if isfield(options, 'ftarget')
    if ~isrealscalar(options.ftarget) || isnan(options.ftarget)
        wid = sprintf('%s:InvalidFtarget', invoker);
        wmsg = sprintf('%s: invalid ftarget; it should be a real number; it is set to %f.', invoker, ftarget);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.ftarget has not got a valid value yet
    options.ftarget = ftarget;
end
options.ftarget = double(options.ftarget);

% Validate options.ctol
validated = false;
if isfield(options, 'ctol')
    if ~isrealscalar(options.ctol) || options.ctol < 0 || isnan(options.ctol)
        wid = sprintf('%s:InvalidCtol', invoker);
        wmsg = sprintf('%s: invalid ctol; it should be a nonnegative number; it is set to %f.', invoker, ctol);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated
    options.ctol = ctol;
end

% Validate options.classical
validated = false;
if isfield(options, 'classical')
    if ~islogicalscalar(options.classical)
        wid = sprintf('%s:InvalidClassicalFlag', invoker);
        wmsg = sprintf('%s: invalid classical flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(classical));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.classical has not got a valid value yet
    options.classical = classical;
end
options.classical = logical(options.classical);
if options.classical
    wid = sprintf('%s:Classical', invoker);
    wmsg = sprintf('%s: in classical mode, which is recommended only for research purpose; set options.classical=false to disable classical mode.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Validate options.fortran
validated = false;
if isfield(options, 'fortran')
    if ~islogicalscalar(options.fortran)
        wid = sprintf('%s:InvalidFortranFlag', invoker);
        wmsg = sprintf('%s: invalid fortran flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(fortran));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif ~options.fortran && options.classical
        wid = sprintf('%s:FortranContradictClassical', invoker);
        wmsg = sprintf('%s: fortran = false but classical = true; fortran is reset to true.', invoker);
        options.fortran = false;
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
        validated = true;
    else
        validated = true;
    end
end
if ~validated % options.fortran has not got a validated value yet
    options.fortran = fortran;
end

% Validate options.honour_x0
validated = false;
if isfield(options, 'honour_x0')
    if ~islogicalscalar(options.honour_x0)
        wid = sprintf('%s:InvalidHonourX0Flag', invoker);
        wmsg = sprintf('%s: invalid honour_x0 flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(honour_x0));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.honour_x0 has not got a valid value yet
    options.honour_x0 = honour_x0;
end
options.honour_x0 = logical(options.honour_x0);

% Validate options.quiet.
validated = false;
% Record user's instruction in the following value; needed for iprint
user_says_quiet = false;
if isfield(options, 'quiet')
    if ~islogicalscalar(options.quiet)
        wid = sprintf('%s:InvalidQuietFlag', invoker);
        wmsg = sprintf('%s: invalid quiet flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(quiet));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
        user_says_quiet = options.quiet;
    end
end
if ~validated % options.quiet has not got a valid value yet
    options.quiet = quiet;
end
options.quiet = logical(options.quiet);

% Validate options.iprint.
validated = false;
if isfield(options, 'iprint')
    if ~isintegerscalar(options.iprint) || (options.iprint ~= 0 && abs(options.iprint) ~= 1 && abs(options.iprint) ~=2 && abs(options.iprint) ~= 3)
        wid = sprintf('%s:InvalidIprint', invoker);
        wmsg = sprintf('%s: invalid iprint; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to %d.', invoker, options.iprint);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif options.iprint ~= 0 && options.classical
        % iprint ~= 0 is not supported in the classical mode.
        wid = sprintf('%s:IprintContradictClassical', invoker);
        wmsg = sprintf('%s: iprint = %d is not supported by the classical mode; it is reset to 0.', invoker, options.iprint);
        options.iprint = 0;
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
        validated = true;
    elseif options.iprint > 0 && user_says_quiet
        % The user says "quiet!" but still asks for information. Let's compromise.
        wid = sprintf('%s:IprintContradictQuiet', invoker);
        if options.classical
            wmsg = sprintf('%s: iprint = %d but quiet = true; iprint is reset to 0.', invoker, options.iprint);
            options.iprint = 0;
            warning(wid, '%s', wmsg);
            warnings = [warnings, wmsg];
            validated = true;
        else
            % In the non-classical mode, we set options.iprint = -options.iprint,
            % meaning that the output will not be displayed on standard output
            % but recorded in a text file SOLVER_output.txt, where SOLVER will
            % be replaced by the solver name. We do not raise a warning since it
            % is explained in the help information and since the user says quiet!
            options.iprint = -options.iprint;
            validated = true;
        end
    elseif options.iprint > 0 && options.fortran
        % iprint > 0 is not supported when calling the Fortran code.
        % This is because of I/O confliction between Fortran and MATLAB.
        wid = sprintf('%s:IprintContradictFortran', invoker);
        wmsg = sprintf('%s: iprint = %d but fortran = true; iprint is reset to %d and the output will be recorded in a .txt file.', invoker, options.iprint, -options.iprint);
        options.iprint = -options.iprint;
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
        validated = true;
    else
        validated = true;
    end
end
if ~validated % options.iprint has not got a valid value yet
    if user_says_quiet
        % The user says "quiet!". Set options.iprint = 0 regarless of the default iprint.
        options.iprint = 0;
    else
        options.iprint = iprint;
    end
end

% Validate options.debug
validated = false;
if isfield(options, 'debug')
    if ~islogicalscalar(options.debug)
        wid = sprintf('%s:InvalidDebugflag', invoker);
        wmsg = sprintf('%s: invalid debug flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(debugflag));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.debug has not got a valid value yet
    options.debug = debugflag;
end
options.debug = logical(options.debug);
if options.debug
    wid = sprintf('%s:Debug', invoker);
    wmsg = sprintf('%s: in debug mode; set options.debug=false to disable debug.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
%    if options.quiet
%        options.quiet = false;
%        wid = sprintf('%s:Debug', invoker);
%        wmsg = sprintf('%s: options.quiet is set to false because options.debug=true.', invoker);
%        warning(wid, '%s', wmsg);
%        warnings = [warnings, wmsg];
%    end
end

% Validate options.chkfunval
validated = false;
if isfield(options, 'chkfunval')
    if ~islogicalscalar(options.chkfunval)
        wid = sprintf('%s:InvalidChkfunval', invoker);
        wmsg = sprintf('%s: invalid chkfunval flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(chkfunval&&options.debug));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif logical(options.chkfunval) && ~options.debug
        wid = sprintf('%s:InvalidChkfunval', invoker);
        wmsg = sprintf('%s: chkfunval=true but debug=false; chkfunval is set to false; set both flags to true to check function values.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.chkfunval has not got a valid value yet
    options.chkfunval = logical(chkfunval) && options.debug;
end
if options.chkfunval
    wid = sprintf('%s:Chkfunval', invoker);
    if strcmp(solver, 'cobylan')
        wmsg = sprintf('%s: checking whether fx=fun(x) and constr=con(x) at exit, which costs an extra function/constraint evaluation; set options.chkfunval=false to disable the check.', invoker);
    else
        wmsg = sprintf('%s: checking whether fx=fun(x) at exit, which costs an extra function evaluation; set options.chkfunval=false to disable the check.', invoker);
    end
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Validate options.maxhist
validated = false;
if isfield(options, 'maxhist')
    if ~isintegerscalar(options.maxhist) || options.maxhist < 0
        wid = sprintf('%s:InvalidMaxhist', invoker);
        wmsg = sprintf('%s: invalid maxhist; it should be a nonnegative integer; it is set to maxfun.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated  % options.maxhist has not got a valid value
    options.maxhist = options.maxfun;  % options.maxfun has been validated
end

% Validate options.output_xhist
validated = false;
if isfield(options, 'output_xhist')
    if ~islogicalscalar(options.output_xhist)
        wid = sprintf('%s:InvalidOutput_xhist', invoker);
        wmsg = sprintf('%s: invalid output_xhist flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(output_xhist));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated
    options.output_xhist = output_xhist;
end

% Validate options.output_nlchist
validated = false;
if isfield(options, 'output_nlchist')
    if ~islogicalscalar(options.output_nlchist)
        wid = sprintf('%s:InvalidOutput_nlchist', invoker);
        wmsg = sprintf('%s: invalid output_nlchist flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(output_nlchist));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated
    options.output_nlchist = output_nlchist;
end

% Validate options.maxfilt
validated = false;
if isfield(options, 'maxfilt')
    if ~isintegerscalar(options.maxfilt) || options.maxfilt < 1
        wid = sprintf('%s:InvalidMaxfilt', invoker);
        wmsg = sprintf('%s: invalid maxfilt; it should be a positive integer; it is set to %d.', invoker, maxfilt);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif options.maxfilt < min_maxfilt
        wid = sprintf('%s:MaxfiltTooSmall', invoker);
        wmsg = sprintf('%s: maxfilt is too small; it should be an integer at least %d; it is set to %d.', invoker, min_maxfilt, maxfilt);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated  % options.maxfilt has not got a valid value
    options.maxfilt = maxfilt;  % options.maxfun has been validated
end

% Validate options.eta1
user_eta1_correct = false;  % Does the user provide a correct eta1? Needed when validating eta2.
validated = false;
if isfield(options, 'eta1')
    if ~isrealscalar(options.eta1) || options.eta1 < 0 || options.eta1 >= 1
        wid = sprintf('%s:InvalidEta1', invoker);
        if isfield(options, 'eta2') && isrealscalar(options.eta2) && options.eta2 > 0 && options.eta2 <= 1
        % The user provides a correct eta2; we define eta1 as follows.
            options.eta1 = max(eps, options.eta2/7);
            wmsg = sprintf('%s: invalid eta1; it should be in the interval [0, 1) and not more than eta2; it is set to %f.', invoker, options.eta1);
            validated = true;
        else
        % The user does not provide a correct eta2; we take the default eta1 hard coded in Powell's code.
            wmsg = sprintf('%s: invalid eta1; it should be in the interval [0, 1) and not more than eta2; it will be set to the default value.', invoker);
        end
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        user_eta1_correct = true;
        validated = true;
    end
end
if ~validated
    options.eta1 = NaN;  % NaN means that Fortran will take the hard-coded default value.
end

% Validate options.eta2
validated = false;
if isfield(options, 'eta2')
    if ~isrealscalar(options.eta2) || (isnan(options.eta1) && options.eta2 < 0) || options.eta2 < options.eta1 || options.eta2 > 1
        wid = sprintf('%s:InvalidEta2', invoker);
        if user_eta1_correct
        % The user provides a correct eta1; we define eta2 as follows.
            options.eta2 = (options.eta1 + 2)/3;
            validated = true;
            wmsg = sprintf('%s: invalid eta2; it should be in the interval [0, 1] and not less than eta1; it is set to %f.', invoker, options.eta2);
        else
        % The user does not provide a correct eta1; we take the default eta2 hard coded in Powell's code.
            wmsg = sprintf('%s: invalid eta2; it should be in the interval [0, 1] and not less than eta1; it will be set to the default value.', invoker);
        end
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated
    options.eta2 = NaN;  % NaN means that Fortran will take the hard-coded default value.
end

% Validate options.gamma1
validated = false;
if isfield(options, 'gamma1')
    if ~isrealscalar(options.gamma1) || options.gamma1 <= 0 || options.gamma1 >= 1
        wid = sprintf('%s:InvalidGamma1', invoker);
        wmsg = sprintf('%s: invalid gamma1; it should be in the interval (0, 1); it will be set to the default value.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated
    options.gamma1 = NaN;  % NaN means that Fortran will take the hard-coded default value.
end

% Validate options.gamma2
validated = false;
if isfield(options, 'gamma2')
    if ~isrealscalar(options.gamma2) || options.gamma2 < 1 || options.gamma2 >= inf
        wid = sprintf('%s:InvalidGamma2', invoker);
        wmsg = sprintf('%s: invalid gamma2; it should be a real number not less than 1; it will be set to the default value.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated
    options.gamma2 = NaN;  % NaN means that Fortran will take the hard-coded default value.
end

% pre_options finished
return

%%%%%%%%%%%%%%%%%%%%%% Function for scaling the problem %%%%%%%%%%%%%%%%
function [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, scaling_factor, shift, substantially_scaled, warnings] = scale_problem(invoker, fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, warnings)
% x_before_scaling = scaling_factor.*x_after_scaling + shift

% Question: What about scaling according to the magnitude of x0, lb, ub,
% x0-lb, ub-x0?
% This can be useful if lb and ub reflect the nature of the problem
% well, and x0 is a reasonable approximation to the optimal solution.
% Otherwise, it may be a bad idea.

callstack = dbstack;
funname =callstack(1).name; % Name of the current function

substantially_scaled_threshold = 2;
% We consider the problem substantially scaled_threshold if
% max([1; scaling_factor])/min([1; scaling_factor]) > substantially_scaled_threshold

% Zaikun 2020-05-24: we change the scaling strategy; do not scale the problem
% unless all variables have both lower and upper bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lenx0 = length(x0);
%scaling_factor = ones(lenx0, 1);
%shift = zeros(lenx0, 1);
%index_lub = (lb > -inf) & (ub < inf); % Variables with lower and upper bounds
%scaling_factor(index_lub) = (ub(index_lub) - lb(index_lub))/2;
%shift(index_lub) = (lb(index_lub) + ub(index_lub))/2;
%shift(~index_lub) = x0(~index_lub); % Shift x0 to 0 unless both lower and upper bounds are present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if max(ub-lb) >= inf  % At least one of [-lb; ub] is infinity
    % Private/unexpcted error
    error(sprintf('%s:InvalidScaling', funname), '%s: UNEXPECTED ERROR: at least one of [-lb; ub] is infinity. Scaling should not be performed.', funname);
end

scaling_factor = (ub - lb)/2;
shift = (lb + ub)/2;

fun = @(x) fun(scaling_factor.*x+shift);
x0 = (x0-shift)./scaling_factor;
if ~isempty(Aineq)
% Aineq*x_before_scaling <= bineq
% <==> Aineq*(scaling_factor.*x_after_scaling+shift) <= bineq
% <==> (Aineq*diag(scaling_factor))*x_after_scaling <= bineq - Aineq*shift
    bineq = bineq - Aineq*shift;
    Aineq = Aineq*diag(scaling_factor);
end
if ~isempty(Aeq)
    beq = beq - Aeq*shift;
    Aeq = Aeq*diag(scaling_factor);
end
if ~isempty(lb)
% lb < x_before_scaling < ub
% <==> lb < scaling_factor.*x_after_scaling + shift < ub
% <==> (lb-shift)./scaling_factor < x_after_scaling < (ub-shift)./scaling_facor
    lb = (lb-shift)./scaling_factor;
end
if ~isempty(ub)
    ub = (ub-shift)./scaling_factor;
end
if ~isempty(nonlcon)
    nonlcon = @(x) nonlcon(scaling_factor.*x+shift);
end

% Zaikun 2020-05-25: We do not warn about scaling anymore. Scaling works
% well in several real problems.
%if any(scaling_factor ~= 1)
%    wid = sprintf('%s:ProblemScaled', invoker);
%    wmsg = sprintf('%s: problem scaled according to bound constraints; do this only if the bounds reflect the scaling of variables; if not, set options.scale=false to disable scaling.', invoker);
%    warning(wid, '%s', wmsg);
%    warnings = [warnings, wmsg];
%end

substantially_scaled = false;
%if (max([scaling_factor; 1./scaling_factor]) > substantially_scaled_threshold)
if max([1; scaling_factor])/min([1; scaling_factor]) > substantially_scaled_threshold
    substantially_scaled = true;
end

if min(scaling_factor) < eps
    % Private/unexpcted error
    error(sprintf('%s:InvalidScaling', funname), '%s: UNEXPECTED ERROR: invalid scaling factor returned when called by %s.', funname, invoker);
end
return

%%%%%%%%%%%%%%%%%%%%%%%% Function for selecting solver %%%%%%%%%%%%%%%%%%%%
function [options, warnings] = select_solver(invoker, options, probinfo, warnings)

invoker_list = {'pdfon'};
% Only pdfo needs select_solver. We may have other invokers in the future!
solver_list = {'uobyqan', 'newuoan', 'bobyqan', 'lincoan', 'cobylan'};
% We may add other solvers in the future!
% Note that pdfo is not a possible solver here!
callstack = dbstack;
funname =callstack(1).name; % Name of the current function

if ~ismember(invoker, invoker_list)
    % Private/unexpected error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s serves only %s.', funname, funname, mystrjoin(invoker_list, ', '));
end
% After pre_options, options.solver is either a member of solver_list
% or '' (i.e., an empty char array), the second signifying the solver
% is yet to decide.
% 1. If options.solver is in solver_list, we check whether it can solve the
% problem. If yes, we set solver=options.solver; otherwise, we warn about
% 'invalid solver' and select a solver.
% 2. If options.solver is '', we do not complain but select a solver. We
% should not complain because either the user does not specify a solver, which
% is perfectly fine, or an unknown solver was specified, which has already
% invoked a warning in pre_options.

solver = options.solver;
ptype = probinfo.refined_type;
n = probinfo.refined_dim;

% Is the user-defined options.solver correct?
solver_correct = ~isempty(solver) && prob_solv_match(ptype, solver);

if ~solver_correct
    if ~isempty(solver) % Do not complain if options.solver is empty.
        wid = sprintf('%s:InvalidSolver', invoker);
        wmsg = sprintf('%s: %s cannot solve a %s problem; %s will select a solver automatically.', invoker, solver, strrep(ptype, '-', ' '), invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    end
    switch ptype
    case 'unconstrained'
        if (n >= 2 && n <= 8 && options.maxfun >= (n+1)*(n+2)/2 + 1)
            solver = 'uobyqan';
        elseif (options.maxfun <= n+2) % After prepdfo, options.maxfun>=n+2 is ensured. Thus options.maxfun<=n+2 indeed means options.maxfun=n+2
            solver = 'cobylan';
        else
            solver = 'newuoan';  % options.npt will be set later
            % Interestingly, we note in our test that lincoan outperformed
            % newuoan on unconstrained CUTEst problems when the dimension
            % was not large (i.e., <=50) or the precision requirement
            % was not high (i.e., >=1e-5). Therefore, it is worthwhile to
            % try lincoan when an unconstrained problem is given.
            % Nevertheless, for the moment, we set the default solver
            % for unconstrained problems to be newuoan.
        end
    case 'bound-constrained'
        if (options.maxfun <= n+2)
            solver = 'cobylan';
        else
            solver = 'bobyqan';  % options.npt will be set later
        end
    case 'linearly-constrained'
        if (options.maxfun <= n+2)
            solver = 'cobylan';
        else
            solver = 'lincoan';  % options.npt will be set later
        end
    case 'nonlinearly-constrained'
        solver = 'cobylan';
    otherwise
        % Private/unexpected error
        error(sprintf('%s:InvalidProbType', funname), '%s: UNEXPECTED ERROR: unknown problem type ''%s'' received.', funname, ptype);
    end
end

% Revise options.npt according to the selected solver
% Note that pre_options has set options.npt to either a positive integer or NaN.
if ismember(solver, {'newuoan', 'bobyqan', 'lincoan'}) && (isnan(options.npt) || options.npt < n+2 || options.npt > min((n+1)*(n+2)/2, options.maxfun-1))
    options.npt = min(2*n+1, options.maxfun - 1);
    if ismember('npt', probinfo.user_options_fields)
        wid = sprintf('%s:InvalidNpt', invoker);
        wmsg = sprintf('%s: npt is set to %d according to the selected solver %s, which requires n+2 <= npt <= (n+1)*(n+2)/2.', invoker, options.npt, solver);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    end
end

% Revise options.rhobeg and options.rhoend according to the selected solver.
% For the moment, only BOBYQA needs such a revision.
if strcmp(solver, 'bobyqan') && options.rhobeg > min(probinfo.refined_data.ub-probinfo.refined_data.lb)/2
    options.rhobeg = max(eps, min(probinfo.refined_data.ub-probinfo.refined_data.lb)/4);
    options.rhoend = max(eps, min(0.1*options.rhobeg, options.rhoend));
    if ismember('rhobeg', probinfo.user_options_fields) || ismember('rhoend', probinfo.user_options_fields)
        wid = sprintf('%s:InvalidRhobeg', invoker);
        wmsg = sprintf('%s: rhobeg is set to %f and rhoend to %f acccording to the selected solver bobyqan, which requires rhoend <= rhobeg <= min(ub-lb)/2.', invoker, options.rhobeg, options.rhoend);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    end
end

if ~ismember(solver, solver_list) || ~prob_solv_match(ptype, solver)
    % Private/unexpected error
    error(sprintf('%s:InvalidSolver', funname), '%s: UNEXPECTED ERROR: invalid solver ''%s'' selected.', funname, solver);
end
options.solver = solver; % Record the solver in options.solver
return

%%%%%%%%%%%%%%%%%%%%%%% Function for checking problem type %%%%%%%%%%%%%%
function ptype = problem_type(Aineq, Aeq, lb, ub, nonlcon)
callstack = dbstack;
funname = callstack(1).name; % Name of the current function

ptype_list = {'unconstrained', 'bound-constrained', 'linearly-constrained', 'nonlinearly-constrained'};

if ~isempty(nonlcon)
    ptype = 'nonlinearly-constrained';
elseif ~isempty(Aineq) || ~isempty(Aeq)
    ptype = 'linearly-constrained';
elseif (~isempty(lb) && max(lb) > -inf) || (~isempty(ub) && min(ub) < inf)
    ptype = 'bound-constrained';
else
    ptype = 'unconstrained';
end

if ~ismember(ptype, ptype_list)
    % Private/unexpected error
    error(sprintf('%s:InvalidProbType', funname), ...
        '%s: UNEXPECTED ERROR: unknown problem type ''%s'' returned.', funname, ptype);
end
return

%%%%%%% Function for checking whether problem type matches solver  %%%%%%
function match = prob_solv_match(ptype, solver)
callstack = dbstack;
funname = callstack(1).name; % Name of the current function

solver_list = {'uobyqan', 'newuoan', 'bobyqan', 'lincoan', 'cobylan', 'pdfon'};
% Note: pdfo is also a possible solver when prob_solv_match is called in
% prepdfo to check whether the invoker can handle the problem.
if ~ismember(solver, solver_list)
    % Private/unexpected error
    error(sprintf('%s:InvalidSolver', funname), ...
    '%s: UNEXPECTED ERROR: unknown solver ''%s'' received.', funname, solver);
end

match = true;
switch ptype
case 'unconstrained'
    match = true;
    % Essentially do nothing. DO NOT remove this case. Otherwise, the
    % case would be included in 'otherwise', which is not correct.
case 'bound-constrained'
    if any(strcmp(solver, {'uobyqan', 'newuoan'}))
        match = false;
    end
case 'linearly-constrained'
    if any(strcmp(solver, {'uobyqan', 'newuoan', 'bobyqan'}))
        match = false;
    end
case 'nonlinearly-constrained'
    if any(strcmp(solver, {'uobyqan', 'newuoan', 'bobyqan', 'lincoan'}))
        match = false;
    end
otherwise
    % Private/unexpected error
    error(sprintf('%s:InvalidProbType', funname), '%s: UNEXPECTED ERROR: unknown problem type ''%s'' received.', funname, ptype);
end
return

%%%%%%%%%%%%%%% Function for calculating constraint violation %%%%%%%%%%
function [constrviolation, nlcineq, nlceq]= constrv(x, Aineq, bineq, Aeq, beq, lb, ub, nonlcon)
% CONSTRV calculates the absolute constraint violation at x.
rineq = [];
req = [];
nlcineq = [];
nlceq = [];
if ~isempty(Aineq)
    rineq = Aineq*x-bineq;
end
if ~isempty(Aeq)
    req = Aeq*x-beq;
end
if ~isempty(nonlcon)
    [nlcineq, nlceq] = nonlcon(x);
end
constrviolation = max([0; rineq; abs(req); lb-x; x-ub; nlcineq; abs(nlceq)], [], 'includenan');
% max(X, [], 'includenan') will return NaN if X contains NaN and the
% maximum of X otherwise.
return

%%%%%% Function for revising x0 or rhobeg when the solver is BOBYQA %%%%
function [x0, options, warnings] = pre_rhobeg_x0(invoker, x0, lb, ub, user_options_fields, options, warnings)
% The Fortran code of BOBYQA will revise x0 so that the distance between x0
% and the inactive bounds is at least rhobeg. We do the revision here in
% order to raise a warning when such a revision occurs. The revision scheme
% is slightly different from the one by Powell in his Fortran code, which sets
% x0 (lb < x0 < lb + rhobeg) = lb + rhobeg
% x0 (ub > x0 > ub - rhobeg) = ub - rhobeg
% Note that lb <= x0 <= ub and rhobeg <= (ub-lb)/2 after pre_options and project.
callstack = dbstack;
funname = callstack(1).name; % Name of the current function

solver_list = {'bobyqan'}; % Only BOBYQA needs pre_rhobeg_x0. May have others in the future.

if ~ismember(lower(options.solver), solver_list)
    % Private/unexpcted error
    error(sprintf('%s:InvalidSolver', funname), '%s: UNEXPECTED ERROR: %s serves only %s.', funname, funname, mystrjoin(solver_list, ', '));
end

if isfield(options, 'honour_x0') && options.honour_x0  % In this case, we respect the user-defiend x0 and revise rhobeg
    rhobeg_old = options.rhobeg;
    lbx = (lb > -inf & x0 - lb <= eps*max(abs(lb), 1));  % x0 essentially equals lb
    ubx = (ub < inf & x0 - ub >= - eps*max(abs(ub), 1));  % x0 essentially equals ub
    options.rhobeg = max(eps, min([options.rhobeg; x0(~lbx) - lb(~lbx); ub(~ubx) - x0(~ubx)]));
    x0(lbx) = lb(lbx);
    x0(ubx) = ub(ubx);
    if rhobeg_old - options.rhobeg > eps*max(1, rhobeg_old)
        options.rhoend = max(eps, min(0.1*options.rhobeg, options.rhoend));  % We do not revise rhoend unless rhobeg is revised
        if ismember('rhobeg', user_options_fields) || ismember('rhoend', user_options_fields)
            wid = sprintf('%s:ReviseRhobeg', invoker);
            wmsg = sprintf('%s: rhobeg is revised to %f and rhoend to %f so that the distance between x0 and the inactive bounds is at least rhobeg.', invoker, options.rhobeg, options.rhoend);
            warning(wid, '%s', wmsg);
            warnings = [warnings, wmsg];
        end
    end
else
    x0_old = x0;
    lbx = (x0 <= lb + 0.5*options.rhobeg);
    lbx_plus = (x0 > lb + 0.5*options.rhobeg) & (x0 < lb + options.rhobeg);
    ubx_minus = (x0 < ub - 0.5*options.rhobeg) & (x0 > ub - options.rhobeg);
    ubx = (x0 >= ub - 0.5*options.rhobeg);
    x0(lbx) = lb(lbx);
    x0(lbx_plus) = lb(lbx_plus) + options.rhobeg;
    x0(ubx_minus) = ub(ubx_minus) - options.rhobeg;
    x0(ubx) = ub(ubx);
    if norm(x0_old-x0) > eps*max(1, norm(x0_old))
        wid = sprintf('%s:ReviseX0', invoker);
        wmsg = sprintf('%s: x0 is revised so that the distance between x0 and the inactive bounds is at least rhobeg; set options.honour_x0=true if you prefer to keep x0.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    end
end
return
