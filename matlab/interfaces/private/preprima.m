function [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options, probinfo] = preprima(fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options)
%PREPRIMA preprocesses the input to prima and its solvers.
%
%   ***********************************************************************
%   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%
%   Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: private (not supposed to be called by users)
%
% Remarks
% 1. Input/output names: MATLAB allows to use the same name for inputs and outputs.
% 2. invoker: invoker is the function that calls preprima.
%
% TODO: None
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preprima starts

warnings = {}; % A cell that records all the warnings, will be recorded in probinfo

% Who is calling this function? Is it a correct invoker?
invoker_list = [all_solvers(), 'prima'];
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if (length(callstack) == 1 || ~ismember(callstack(2).name, invoker_list))
    % Private/unexpected error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s should only be called by %s.', funname, funname, strjoin(invoker_list, ', '));
else
    invoker = callstack(2).name; % Name of the function who calls this function
end

if (nargin ~= 1) && (nargin ~= 10)
    % Private/unexpected error
    error(sprintf('%s:InvalidInput', funname), '%s: UNEXPECTED ERROR: 1 or 10 inputs.', funname);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If invoker is a solver called by prima, then preprima should have been called in prima.
if (length(callstack) >= 3) && strcmp(callstack(3).name, 'prima')
    if nargin ~= 10 % There should be 10 input arguments
        % Private/unexpected error
        error(sprintf('%s:InvalidInput', funname), ...
        '%s: UNEXPECTED ERROR: %d inputs received; this should not happen as preprima has been called once in prima.', funname, nargin);
    end
    % In this case, we set probinfo to empty.
    probinfo = [];
    return % Return because preprima has already been called in prima.
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
% 10. trivial_leq: a true/false vector indicating which linear equality
%     constraints are trivial (up to naive tests)
% 11. infeasible: whether the problem is infeasible (up to naive tests)
% 12. scaled: whether the problem is scaled
% 13. scaling_factor: vector of scaling factors
% 14. shift: vector of shifts
% 15. reduced: whether the problem is reduced (due to fixed variables)
% 16. raw_type: problem type before reduction
% 17. raw_dim: problem dimension before reduction
% 18. refined_type: problem type after reduction
% 19. refined_dim: problem dimension after reduction
% 20. feasibility_problem: whether the problem is a feasibility problem
% 21. user_options_fields: the fields in the user-specified options
% 22. options: (refined) options for calling the solvers
% 23. warnings: warnings during the preprocessing/validation
% 24. boundmax: the large allowed absolute value of bound constraints
% 25. funcmax: the largest allowed value of the objective function
% 26. constrmax: the largest allowed absolute value of the constraints
% 27. x0_is_row: whether x0 is a row
probinfo = struct();

% Save the raw data (date before validation/preprocessing) in probinfo.
% The raw data can be useful when debugging. At the end of preprima, if
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

% Decide the precision ('single', 'double', or 'quadruple') of the real calculation within the
% Fortran solvers. This is needed ONLY by `boundmax`, `funcmax`, and `constrmax` defined below by
% calling `getmax`. These three numbers will be used in `pre_x0`, `pre_fun`, and `pre_nonlcon`
% respectively in the sequel. Note the following.
% 1. `precision` takes effect only if Fortran solvers are called (i.e., when options.fortran = true).
% 2. `precision` is passed only to `getmax`, which defines huge values (e.g., `boundmax`, `funcmax`).
precision = 'double';
% Since `options` is not validated yet, validations are needed before inquiring options.precision.
if isa(options, 'struct') && isfield(options, 'precision') && ischarstr(options.precision) && ...
        ismember(lower(options.precision), all_precisions())
    precision = lower(options.precision);
end
probinfo.boundmax = getmax('bound', precision);
probinfo.funcmax = getmax('function', precision);
probinfo.constrmax = getmax('constraint', precision);

% Validate and preprocess x0
[x0, warnings, x0_is_row] = pre_x0(invoker, x0, precision, warnings);
lenx0 = length(x0); % Within this file, for clarity, we denote length(x0) by lenx0 instead of n
probinfo.x0_is_row = x0_is_row;

% Validate and preprocess the bound constraints
% In addition, get the indices of infeasible bounds and 'fixed variables'
% such that ub-lb < 2*eps (if any) and save the information in probinfo.
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
% min(ub) and max(lb) are evaluated) and before preprocessing the
% linear/nonlinear constraints (because these constraints will be
% reduced during the preprocessing). Note that Aineq, Aeq, and nonlcon will
% not be "evaluated" in problem_type, so there is no worry about the
% validity of them.
probinfo.raw_type = problem_type(Aineq, Aeq, lb, ub, nonlcon);

% Raise a warning if x0 is a row and the problem has linear constraints (before the reduction),
% as the formulation of linear constraints (Aineq*x <= bineq, Aeq*x = beq) assumes that the decision
% variable x is a column. If x0 is a row, pre_x0 transposes it, and pre_fun, pre_nonlcon redefines
% fun and nonlcon, so that the solvers only need to deal with x as a column.
if x0_is_row && ~(isempty(Aineq) && isempty(Aeq))
    wid = sprintf('%s:X0IsRow', invoker);
    wmsg = sprintf('%s: x0 is a row, but a column is expected; it is transposed.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Validate and preprocess the linear constraints
% 1. The constraints will be reduced if some but not all variables are
%    fixed by the bound constraints. See pre_lcon for why we do not
%    reduce the problem when all variables are fixed.
% 2. The 'trivial constraints' will be excluded (if any).
% 3. In addition, get the indices of infeasible and trivial constraints (if any)
%    and save the information in probinfo.
% 4. fixedx_value is revised to the corresponding values of x0 if infeasible linear constraints are detected
[Aineq, bineq, Aeq, beq, infeasible_lineq, trivial_lineq, infeasible_leq, trivial_leq, fixedx_value, warnings] = pre_lcon(invoker, x0, Aineq, bineq, Aeq, beq, lenx0, fixedx, fixedx_value, warnings);
probinfo.fixedx_value = fixedx_value; % Possibly revised value of the fixed x entries
probinfo.infeasible_lineq = infeasible_lineq; % A vector of true/false
probinfo.trivial_lineq = trivial_lineq; % A vector of true/false
probinfo.infeasible_leq = infeasible_leq; % A vector of true/false
probinfo.trivial_leq = trivial_leq; % A vector of true/false

% Validate and preprocess fun
% 1. The objective function will be reduced if some but not all variables are
%    fixed by the bound constraints. See pre_lcon for why we do not reduce the
%    problem when all variables are fixed.
% 2. This should be done after preprocessing the bound and linear constraints, which defines
%    and may revise fixedx_value.
[fun, probinfo.feasibility_problem, warnings] = pre_fun(invoker, fun, fixedx, fixedx_value, probinfo.funcmax, x0_is_row, warnings);

% Validate and preprocess the nonlinear constraints
% 1. The constraints will be reduced if some but not all variables are fixed by the bound
%    constraints. See pre_lcon for why we do not reduce the problem when all variables
%    are fixed.
% 2. This should be done after preprocessing the bound and linear constraints, which defines
%    and may revise fixedx_value.
% 3. This should be done before evaluating probinfo.constrv_x0 or probinfo.constrv_fixedx.
nonlcon = pre_nonlcon(invoker, nonlcon, fixedx, fixedx_value, probinfo.constrmax, x0_is_row);

% Reduce x0, lb, and ub if some but not all variables are fixed by
% the bound constraints. See pre_lcon for why we do not reduce the
% problem when all variables are fixed.
probinfo.raw_dim = lenx0; % Problem dimension before reduction
if any(fixedx) && any(~fixedx)
    freex = ~fixedx; % A vector of true/false indicating whether the variable is free or not
    x0 = x0(freex); % x0 after reduction
    lenx0 = length(x0);
    lb = lb(freex); % lb after reduction
    ub = ub(freex); % ub after reduction
end
probinfo.refined_type = problem_type(Aineq, Aeq, lb, ub, nonlcon); % Problem type after reduction
probinfo.refined_dim = length(x0); % Problem dimension after reduction
probinfo.reduced = any(fixedx) && any(~fixedx); % Whether the problem has been reduced

% After the preprocessing, the problem may turn out infeasible
probinfo.infeasible = any([probinfo.infeasible_lineq; probinfo.infeasible_leq; probinfo.infeasible_bound]);
if probinfo.infeasible  % The problem turns out infeasible
    [probinfo.constrv_x0, probinfo.nlcineq_x0, probinfo.nlceq_x0] = get_constrv(x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon);
    % The constraint violation calculated by constrv does not include
    % the violation of x0 for the bounds corresponding to fixedx; the
    % corresponding values of x0 are in fixedx_value, while the values
    % of the bounds (lb and ub are the same up to eps) are in
    % fixedx_value_save. Thus the violation is abs(fixedx_value-fixedx_value_save).
    rbounds = abs(fixedx_value - fixedx_value_save);
    rbounds = rbounds(fixedx_value ~= fixedx_value_save);  % Prevent NaN in case both are +/-Inf
    probinfo.constrv_x0 = max([probinfo.constrv_x0; rbounds], [], 'includenan');
end

% After the preprocessing, x may turn out fixed by the bounds
probinfo.nofreex = all(fixedx);
if probinfo.nofreex  % x turns out fixed by the bound constraints
    [probinfo.constrv_fixedx, probinfo.nlcineq_fixedx, probinfo.nlceq_fixedx] = get_constrv(probinfo.fixedx_value, Aineq, bineq, Aeq, beq, lb, ub, nonlcon);
end

% Can the invoker handle the given problem?
% This should be done after the problem type has bee 'refined'.
if ~prob_solv_match(probinfo.refined_type, invoker)
    if strcmp(invoker, 'prima') || (nargin ~= 1)
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

% Revise x0 for bound and linearly constrained problems.
% This is necessary for LINCOA, which accepts only feasible x0.
% Should we do this even if there are nonlinear constraints?
% For now, we do not, because doing so may dramatically increase the
% infeasibility of x0 with respect to the nonlinear constraints.
if ismember(probinfo.refined_type, {'bound-constrained', 'linearly-constrained'}) && ~probinfo.nofreex && ~probinfo.infeasible
    % Another possibility for bound-constrained problems:
    % xind = (x0 < lb) | (x0 > ub);
    % x0(xind) = (lb(xind) + ub(xind))/2;
    x0_new = project(Aineq, bineq, Aeq, beq, lb, ub, x0);
    if get_cstrv(x0_new, Aineq, bineq, Aeq, beq, lb, ub) < get_cstrv(x0, Aineq, bineq, Aeq, beq, lb, ub)
        if norm(x0_new-x0) > eps*max(1, norm(x0)) && ~probinfo.feasibility_problem
            % No warning about revising x0 if the problem is a linear feasibility problem.
            % Note that the linearity is guaranteed by THE OUTER IF.
            wid = sprintf('%s:ReviseX0', invoker);
            wmsg = sprintf('%s: x0 is revised to satisfy the constraints better.', invoker);
            warning(wid, '%s', wmsg);
            warnings = [warnings, wmsg];
        end
        x0 = x0_new;
    end
    if get_cstrv(x0, Aineq, bineq, Aeq, beq, lb, ub) > 1.0e-10*max(abs([1; bineq; beq; x0]))
        wid = sprintf('%s:InfeasibleX0', invoker);
        wmsg = sprintf('%s: preprocessing code did not find a feasible x0; problem is likely infeasible.', invoker);
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
% and probinfo.refined_data.ub will be used for defining rhobeg if bobyqa is selected.
probinfo.refined_data = struct('objective', fun, 'x0', x0, 'Aineq', Aineq, 'bineq', bineq, ...
    'Aeq', Aeq, 'beq', beq, 'lb', lb, 'ub', ub, 'nonlcon', nonlcon);

% Select a solver if invoker = 'prima'; record the solver in options.solver.
% Some options will be revised accordingly, including npt, rhobeg, rhoend.
% Of course, if the user-defined options.solver is valid, we accept it.
if strcmp(invoker, 'prima')
    [options, warnings] = select_solver(invoker, options, probinfo, warnings);
end

% If options.fortran is true, check whether the Fortran MEX function is available. If no, set
% options.fortran to false, and the MATLAB version of the solver will be called (we are assuming
% that the MATLAB version is available).
if options.fortran
    [options, warnings] = is_fortran_available(invoker, options, warnings);
end

if strcmpi(options.solver, 'bobyqa') && ~probinfo.nofreex && ~probinfo.infeasible && ~probinfo.feasibility_problem
% BOBYQA will revise x0 so that the distance between x0 and the inactive bounds
% is at least rhobeg. We do it here in order to raise a warning when such a
% revision occurs. After this, BOBYQA will not revise x0 again. If options.honour_x0
% is true, then we keep x0 unchanged and revise rhobeg if necessary.
% N.B.: If x0 violates the bounds, then it is always revised by `project` to respect the bounds.
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
% When the problem is a linear feasibility problem, PRIMA will return the
% current x0, which has been revised by project. The constraint violation
% at x0 is needed to set the output. Note that there is no nonlinear
% constraint in this case.
    probinfo.constrv_x0 = get_constrv(x0, Aineq, bineq, Aeq, beq, lb, ub, []);
end

probinfo.warnings = warnings; % Record the warnings in probinfo

if ~options.debug % Do not carry the raw data with us unless in debug mode.
    probinfo.raw_data = struct();
    % Set this field to empty instead of remove it, because postprima
    % requires this field to exist.
end

if ~options.debug && ~probinfo.scaled
    % The refined data is used only when the problem is scaled. It can
    % also be useful when debugging.
    probinfo.refined_data = struct();
    % Set this field to empty instead of remove it, because postprima
    % requires this field to exist.
end

% preprima ends
return

%%%%%%%%%%%%%%%%%%%%%%%% Function for problem decoding %%%%%%%%%%%%%%%%%
function [fun, x0, Aineq, bineq, Aeq, beq, lb, ub, nonlcon, options, warnings] = decode_problem(invoker, problem, warnings)
% Read the fields of the 'problem' structure but do not validate them.
% The decoded problem will be sent to the preprima function for validation.
% NOTE: We treat field names case-sensitively.

% Possible invokers
invoker_list = [all_solvers(), 'prima'];

callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if ~ismember(invoker, invoker_list)
    % invoker affects the behavior of this function, so we check invoker
    % again, even though it should have been checked in function preprima
    % Private/unexpected error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s serves only %s.', funname, funname, strjoin(invoker_list, ', '));
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
    '%s: PROBLEM misses the %s field(s).', invoker, strjoin(missing_fields, ', '));
end
x0 = problem.x0;

if isfield(problem, 'objective')
    fun = problem.objective;
else % There is no objective; this is a feasibility problem
    fun = []; % We use [] to signify that fun is not specified. pre_fun will replace [] by by @(x)0
end

% Are there unknown fields?
known_fields = {'objective', 'x0', 'Aineq', 'bineq', 'Aeq', 'beq', 'lb', 'ub', 'nonlcon', 'options', 'solver'};
% 1. When invoker is in {uobyqa, ..., cobyla}, we will not complain that
%    a solver is specified unless invoker ~= solver. See function pre_options.
% 2. When invoker is in {uobyqa, ..., cobyla}, if the problem turns out
%    unsolvable for the invoker, then we will raise an error in preprima.
%    We do not do it here because the problem has not been validated/preprocessed
%    yet. Maybe some constraints are trivial and hence can be removed
%    (e.g., bineq = inf, lb = -inf), which can change the problem type.

unknown_fields = setdiff(problem_fields, known_fields);
problem = rmfield(problem, unknown_fields);  % Remove the unknown fields

if ~isempty(unknown_fields)
    wid = sprintf('%s:UnknownProbField', invoker);
    if length(unknown_fields) == 1
        wmsg = sprintf('%s: problem with an unknown field %s; it is ignored.', invoker, strjoin(unknown_fields, ', '));
    else
        wmsg = sprintf('%s: problem with unknown fields %s; they are ignored.', invoker, strjoin(unknown_fields, ', '));
    end
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Read the fields of problem. They will be validated in function preprima
Aineq = [];
bineq = [];
Aeq = [];
beq = [];
lb = [];
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

%%%%%%%%%%%%%%%%%%%%%%%% Function for x0 preprocessing %%%%%%%%%%%%%%%%%
function [x0, warnings, x0_is_row] = pre_x0(invoker, x0, precision, warnings)
[isrv, lenx0]  = isrealvector(x0);
if ~(isrv && (lenx0 > 0))
    % Public/normal error
    error(sprintf('%s:InvalidX0', invoker), '%s: X0 should be a real vector/scalar.', invoker);
end
if (lenx0 > maxint())
    % Public/normal error
    error(sprintf('%s:ProblemTooLarge', invoker), '%s: The problem is too large; at most %d variables are allowed.', invoker, maxint());
end
x0_is_row = (lenx0 >= 2) && isrealrow(x0);
x0 = double(x0(:));  % Transpose x0 if it is a row; fun and nonlcon will be redefined accordingly.
abnormal_x0 = isnan(x0) | (abs(x0) >= inf);
if any(abnormal_x0)
    x0(isnan(x0)) = 0;
    maxfloat = getmax('real', precision);
    x0(isinf(x0) & x0 > 0) = maxfloat;
    x0(isinf(x0) & x0 < 0) = -maxfloat;
    wid = sprintf('%s:AbnormalX0', invoker);
    wmsg = sprintf('%s: X0 contains NaN or infinite values; NaN is replaced by 0 and Inf by %g.', invoker, maxfloat);
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
    '%s: lb should be a real vector and length(lb) = length(x0) unless lb = [].', invoker);
end
if (lenlb == 0)
    lb = -inf(lenx0,1); % After pre_bcon, length(lb) = length(x0)
end
lb = double(lb(:));
if any(isnan(lb))
    wid = sprintf('%s:NaNInLB', invoker);
    wmsg = sprintf('%s: LB contains NaN; the problem is hence infeasible.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Upper bounds (ub)
[isrvub, lenub] = isrealvector(ub);
if ~(isrvub && (lenub == lenx0 || lenub == 0))
    % Public/normal error
    error(sprintf('%s:InvalidBound', invoker), ...
    '%s: ub should be a real vector and length(ub) = length(x0) unless ub = [].', invoker);
end
if (lenub == 0)
    ub = inf(lenx0,1); % After pre_bcon, length(ub) = length(x0)
end
ub = double(ub(:));
if any(isnan(ub))
    wid = sprintf('%s:NaNInUB', invoker);
    wmsg = sprintf('%s: UB contains NaN; the problem is hence infeasible.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

infeasible_bound = ~(lb <= ub); % A vector of true/false; true if lb or ub is NaN or lb > ub
if any(infeasible_bound)
    fixedx = false(lenx0, 1);
    fixedx_value = [];
else
    fixedx = (ub <= lb + 2*eps);  % Avoid ub - lb in case ub = lb = +/-Inf. Use <= instead of <.
    fixedx_value = (lb(fixedx)+ub(fixedx))/2;
end
return

%%%%%%%%%%%%%%%%% Function for linear constraint preprocessing %%%%%%%%%%
function [Aineq, bineq, Aeq, beq, infeasible_lineq, trivial_lineq, infeasible_leq, trivial_leq, fixedx_value, warnings] = pre_lcon(invoker, x0, Aineq, bineq, Aeq, beq, lenx0, fixedx, fixedx_value, warnings)

freex = ~fixedx; % A vector of true/false indicating whether the variable is free or not

% Preprocess linear inequalities: Aineq*x <= bineq.

% Check whether Aineq and bineq are real matrices and vectors of proper size, respectively.
[isrm, mA, nA] = isrealmatrix(Aineq);
[isrv, lenb] = isrealvector(bineq);  % The same as fmincon, we allow bineq to be a row
% No matter whether x0 or bineq is a row or column, we always require that size(Aineq) is
% [length(bineq), length(x0)] unless Aineq = bineq = []. This is consistent with fmincon.
if ~(isrm && isrv && (mA == lenb) && (nA == lenx0 || nA == 0))
    % Public/normal error
    error(sprintf('%s:InvalidLinIneq', invoker), ...
    '%s: Aineq should be a real matrix, bineq should be a real column, and size(Aineq) = [length(bineq), length(X0)] unless Aineq = bineq = [].', invoker);
end
if (lenb >= 2) && isrealrow(bineq)
    wid = sprintf('%s:BineqIsRow', invoker);
    wmsg = sprintf('%s: bineq is a row, but a column is expected; it is transposed.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end
bineq = double(bineq(:));
Aineq = double(Aineq);

% Warn about inequality constraints containing NaN.
nan_ineq = any(isnan(Aineq), 2) | isnan(bineq);
if any(nan_ineq)
    wid = sprintf('%s:NaNInequality', invoker);
    wmsg = sprintf('%s: Aineq or bineq contains NaN; the problem is hence infeasible.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Warn about inequality constraints containing Inf.
inf_ineq = any(abs(Aineq) >= inf, 2) | (bineq <= -inf);
if any(inf_ineq)
    wid = sprintf('%s:InfInequality', invoker);
    wmsg = sprintf('%s: Aineq contains infinite values or bineq contains -Inf; the problem is considered infeasible.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Reduce the inequality constraints if some but not all variables are fixed.
% 1. This has to be done before detecting the "zero constraints" (i.e., constraints with zero
%    gradients), because nonzero constraints may become zero after reduction.
% 2. We should NOT reduce the problem if all variables are fixed. Otherwise, Aineq would be [], and
%    then bineq will be set to [] in the end. In this way, we lose completely the information in
%    these constraints. Consequently, we cannot evaluate the constraint violation correctly when needed.
lineq_reduced = false; % Whether linear inequality constraints are reduced
if ~isempty(Aineq) && any(fixedx) && any(~fixedx)
    Aineq_fixed = Aineq(:, fixedx); % Aineq_fixed and bineq_save will be used when revising fixedx_value
    bineq_save = bineq;
    bineq = bineq - Aineq_fixed * fixedx_value;
    Aineq = Aineq(:, freex);
    lineq_reduced = true;
end

% Define infeasible_lineq.
if isempty(Aineq)
    infeasible_lineq = [];
else
    Aineq_rownorm1 = sum(abs(Aineq), 2);
    zero_ineq = (Aineq_rownorm1 == 0);
    Aineq_rownorm1(zero_ineq) = 1;
    infeasible_zero_ineq = (bineq < 0 & zero_ineq);
    % bineq has been revised during the reduction; we regard the constraint as infeasible if Inf or
    % NaN arises after the reduction.
    nan_ineq = nan_ineq | isnan(bineq);
    inf_ineq = inf_ineq | (abs(bineq) >= inf);
    infeasible_lineq = (bineq ./ Aineq_rownorm1 <= -inf) | infeasible_zero_ineq | nan_ineq | inf_ineq; % A vector of true/false
end

% Preprocess linear equalities: Aeq*x == beq

% Check whether Aeq and beq are real matrices and vectors of proper size, respectively.
[isrm, mA, nA] = isrealmatrix(Aeq);
[isrv, lenb] = isrealvector(beq);  % The same as fmincon, we allow beq to be a row
% No matter whether x0 or beq is a row or column, we always require that size(Aeq) is
% [length(beq), length(x0)] unless Aeq = beq = []. This is consistent with fmincon.
if ~(isrm && isrv && (mA == lenb) && (nA == lenx0 || nA == 0))
    % Public/normal error
    error(sprintf('%s:InvalidLinEq', invoker), ...
    '%s: Aeq should be a real matrix, beq should be a real column, and size(Aeq) = [length(beq), length(X0)] unless Aeq = beq = [].', invoker);
end
if (lenb >= 2) && isrealrow(beq)
    wid = sprintf('%s:BeqIsRow', invoker);
    wmsg = sprintf('%s: beq is a row, but a column is expected; it is transposed.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end
beq = double(beq(:));
Aeq = double(Aeq);

% Warn about equality constraints containing NaN.
nan_eq = any(isnan(Aeq), 2) | isnan(beq);
if any(nan_eq)
    wid = sprintf('%s:NaNEquality', invoker);
    wmsg = sprintf('%s: Aeq or beq contains NaN; The problem is hence infeasible.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Warn about equality constraints containing Inf.
inf_eq = any(abs(Aeq) >= inf, 2) | (abs(beq) >= inf);
if any(inf_eq)
    wid = sprintf('%s:InfEquality', invoker);
    wmsg = sprintf('%s: Aeq or beq contains infinite values; the problem is considered infeasible.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Reduce the equality constraints if some but not all variables are fixed.
% 1. This has to be done before detecting the "zero constraints" (i.e., constraints with zero
%    gradients), because nonzero constraints may become zero after reduction.
% 2. We should NOT reduce the problem if all variables are fixed. Otherwise, Aeq would be [], and
%    then beq will be set to [] in the end. In this way, we lose completely the information in
%    these constraints. Consequently, we cannot evaluate the constraint violation correctly when needed.
leq_reduced = false; % Whether linear equality constraints are reduced
if ~isempty(Aeq) && any(fixedx) && any(~fixedx)
    Aeq_fixed = Aeq(:, fixedx); % Aeq_fixed and beq_save may be used when revising fixedx_value
    beq_save = beq;
    beq = beq - Aeq_fixed * fixedx_value;
    Aeq = Aeq(:, freex);
    leq_reduced = true;
end

% Define infeasible_leq.
if isempty(Aeq)
    infeasible_leq = [];
else
    Aeq_rownorm1 = sum(abs(Aeq), 2);
    zero_eq = (Aeq_rownorm1 == 0);
    Aeq_rownorm1(zero_eq) = 1;
    infeasible_zero_eq = (beq ~= 0 & zero_eq);
    % beq has been revised during the reduction; we regard the constraint as infeasible if Inf or
    % NaN arises after the reduction.
    nan_eq = nan_eq | isnan(beq);
    inf_eq = inf_eq | (abs(beq) >= inf);
    infeasible_leq = (abs(beq ./ Aeq_rownorm1) >= inf) | infeasible_zero_eq | nan_eq | inf_eq;  % A vector of true/false
end

% Define trivial_lineq and trivial_leq; remove the trivial constraints.
infeasible = (any(infeasible_lineq) || any(infeasible_leq));
if infeasible
    trivial_lineq = false(size(bineq));
    trivial_leq = false(size(beq));
else
    if isempty(Aineq)
        trivial_lineq = [];
    else
        trivial_lineq = ((bineq ./ Aineq_rownorm1 == inf) | (bineq >= inf) | (bineq >= 0 & zero_ineq));
        Aineq = Aineq(~trivial_lineq, :); % Remove the trivial linear inequalities
        bineq = bineq(~trivial_lineq);
    end
    if isempty(Aeq)
        trivial_leq = [];
    else
        trivial_leq = (beq == 0 & zero_eq);
        Aeq = Aeq(~trivial_leq, :); % Remove trivial linear equalities
        beq = beq(~trivial_leq);
    end
end

% If infeasibility is detected, then we will return x0 without further calculations. Thus we need to
% revise fixedx_value to x0(fixedx), and redefine bineq/beq so that they are reduced with
% x(fixedx) = x0(fixedx) (otherwise, the constraint violation cannot be correctly calculated later).
reduced = (lineq_reduced || leq_reduced);
if infeasible && reduced
    fixedx_value = x0(fixedx);
    if lineq_reduced
        bineq = bineq_save - Aineq_fixed * fixedx_value;
    end
    if leq_reduced
        beq = beq_save - Aeq_fixed * fixedx_value;
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

if (max([length(beq), length(bineq)]) > maxint())
    % Public/normal error
    error(sprintf('%s:ProblemTooLarge', invoker), '%s: The problem is too large; at most %d constraints are allowed.', invoker, maxint());
end
return

%%%%%%%%%%%%%%%%%%%%%%%% Function for fun preprocessing %%%%%%%%%%%%%%%%%
function [fun, feasibility_problem, warnings] = pre_fun(invoker, fun, fixedx, fixedx_value, funcmax, x0_is_row, warnings)
if ~(isempty(fun) || ischarstr(fun) || isa(fun, 'function_handle'))
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
elseif ischarstr(fun)
    fun = str2func(fun);
    % Work with function handles instead of function names to avoid using 'feval'
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
% During the calculation, x is always a column vector. We assume that fun expects a row vector if
% x0 is a row.
if x0_is_row
    fun = @(x) evalobj(invoker, fun, x', funcmax);  % See evalobj.m for evalobj
else
    fun = @(x) evalobj(invoker, fun, x, funcmax);  % See evalobj.m for evalobj
end
% Reduce fun if some but not all variables are fixed by the bounds.
% Note that we do not reduce the problem when all variables are fixed. See pre_lcon for the reason.
if any(fixedx) && any(~fixedx)
    fun = @(freex_value) fun(fullx(freex_value, fixedx_value, fixedx));
end
return

%%%%%%%%%%%%%%%%% Function for nonlinear constraint preprocessing %%%%%%%%%%
function nonlcon = pre_nonlcon(invoker, nonlcon, fixedx, fixedx_value, constrmax, x0_is_row)
if ~(isempty(nonlcon) || isa(nonlcon, 'function_handle') || ischarstr(nonlcon))
    % Public/normal error
    error(sprintf('%s:InvalidCon', invoker), ...
    '%s: nonlcon should be a function handle or a function name.', invoker);
end
if isempty(nonlcon)
    nonlcon = []; % We use [] to signify that nonlcon is not specified; its size is 0x0
else
    if ischarstr(nonlcon)
        nonlcon = str2func(nonlcon);
        % work with function handles instead of function names to avoid using 'feval'
    end
    if ~exist('OCTAVE_VERSION', 'builtin')
        % Check whether nonlcon has at least 2 outputs.
        % nargout(fun) = #outputs in the definition of fun.
        % If fun includes varargout in definition, nargout(fun) = -#outputs.
        % Octave does not support nargout for built-in functions (as of 2019-08-16)!
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
            '%s: nonlcon has too few outputs; it should return [cineq, ceq], the constraints being cineq(x) <= 0, ceq(x) = 0.', invoker);
        end
    end
    % During the calculation, x is always a column vector. We assume that nonlcon expects a row
    % vector if x0 is a row.
    if x0_is_row
        nonlcon = @(x) evalcon(invoker, nonlcon, x', constrmax);  % See evalcon.m for evalcon
    else
        nonlcon = @(x) evalcon(invoker, nonlcon, x, constrmax);  % See evalcon.m for evalcon
    end
    % Reduce the nonlcon if some but not all variables are fixed by the bounds
    % Note that we do not reduce the problem when all variables are fixed. See pre_lcon for the reason.
    if any(fixedx) && any(~fixedx)
        nonlcon = @(freex_value) nonlcon(fullx(freex_value, fixedx_value, fixedx));
    end
end
return

%%%%%%%%%%%%%%%%% Function fullx used when reducing the problem %%%%%%%%
function x = fullx(freex_value, fixedx_value, fixedx)
x = NaN(length(freex_value)+length(fixedx_value), 1);
x(~fixedx) = freex_value;
x(fixedx) = fixedx_value;
return

%%%%%%%%%%%%%%%%% Function for option preprocessing %%%%%%%%%%
function [options, user_options_fields, warnings] = pre_options(invoker, options, lenx0, lb, ub, warnings)

% NOTE: We treat field names case-sensitively.

% Possible solvers
solver_list = all_solvers();
% We may add other solvers in the future!
% If a new solver is included, we should do the following.
% 0. Include it into the invoker_list (in this and other functions).
% 1. What options does it expect? Set known_fields accordingly.
% 2. Set default options accordingly.
% 3. Check other functions (especially decode_problem, whose behavior
%    depends on the invoker/solver. See known_fields there).

% Possible invokers
invoker_list = [solver_list, 'prima'];

callstack = dbstack;
funname = callstack(1).name;
% invoker affects the behavior of this function, so we check invoker
% again, even though it should have been checked in function preprima
if ~ismember(invoker, invoker_list)
    % Private/unexpected error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s serves only %s.', funname, funname, strjoin(invoker_list, ', '));
end

% Default values of the options.
% npt = ! LATER ! % The default npt depends on solver and will be set later in this function
maxfun = 500*lenx0;
rhobeg = 1; % The default rhobeg and rhoend will be revised if solver = 'bobyqa'
rhoend = 1e-6;
ftarget = -inf;
ctol = sqrt(eps); % Tolerance for constraint violation; a point with a constraint violation at most ctol is considered feasible
cweight = 1e8;  % The weight of constraint violation in the selection of the returned x
classical = false; % Call the classical Powell code? Classical mode recommended only for research purpose
precision = 'double'; % The precision of the real calculation within the solver
fortran = true; % Call the Fortran code?
scale = false; % Scale the problem according to bounds? Scale only if the bounds reflect well the scale of the problem
scale = (scale && max(ub-lb)<inf); % ! NEVER remove this ! Scale only if all variables are with finite lower and upper bounds
honour_x0 = false; % Respect the user-defined x0? Needed by BOBYQA
iprint = 0;
quiet = true;
debug_flag = false; % Do not use 'debug' as the name, which is a MATLAB function
chkfunval = false;
output_xhist = false; % Output the history of x?
output_nlchist = false; % Output the history of the nonlinear constraints?
min_maxfilt = 200; % The smallest value of maxfilt; if maxfilt is too small, the returned x may not be the best one visited
maxfilt = 10*min_maxfilt; % Length of the filter used for selecting the returned x in constrained problems

if ~(isa(options, 'struct') || isempty(options))
    % Public/normal error
    error(sprintf('%s:InvalidOptions', invoker), '%s: OPTIONS should be a structure.', invoker);
end

% Which fields are specified?
options = rmempty(options); % Remove empty fields
options_fields = fieldnames(options);
% The list of fields in options  will be returned and used elsewhere. We
% save it right now in case we "intelligently" change options_fields
% after this line in future versions.
user_options_fields = options_fields;

% Validate options.solver
% We need to know what is the solver in order to decide which fields
% are 'known' (e.g., expected), and also to set npt, rhobeg, rhoend.
% We do the following:
% 1. If invoker = 'prima':
% 1.1 If no solver is specified or solver = 'prima', we do not complain
% and set options.solver = solver = '', i.e., an empty char array;
% 1.2 Else if solver is not in solver_list, we warn about 'unknown solver'
% and set options.solver = solver = '', i.e., an empty char array;
% 1.3 Else, we set solver = options.solver.
% 2. If invoker is in solver_list:
% 2.1 If options.solver exists but options.solver ~= invoker, we warn
% about 'invalid solver' and set options.solver = solver = invoker;
% 2.2 Else, we do not complain and set options.solver = solver = invoker.
% In this way, options.solver and solver either end up with a member of
% solver_list or ''. The second case is possible only if invoker = 'prima',
% and solver will be selected later.
if isfield(options, 'solver') && ~ischarstr(options.solver)
    options.solver = 'UNKNOWN_SOLVER';
    % We have to change options.solver to a char/string so that we can use strcmpi
    % We do not need to worry about the case where solver is empty, because
    % all the empty fields have been removed from options.
end
if strcmp(invoker, 'prima')
    % We set the default value of solver to '', an empty char array.
    % 1. DO NOT change this default value! It will affect known_fields
    % and select_solver.
    % 2. DO NOT use [], which is an empty double array and may cause some
    % functions (e.g., ismember) to complain about incompatible types.
    solver = '';
    if isfield(options, 'solver')
        if any(strcmpi(options.solver, solver_list))
            solver = lower(options.solver);
        elseif ~strcmpi(options.solver, 'prima')
        % We should not complain about 'unknown solver' if invoker = options.solver = 'prima'
            wid = sprintf('%s:UnknownSolver', invoker);
            wmsg = sprintf('%s: unknown solver specified; %s will select one automatically.', invoker, invoker);
            warning(wid, '%s', wmsg);
            warnings = [warnings, wmsg];
        end
    end
else % invoker is in {'uobyqa', ..., 'cobyla'}
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
options.solver = char(solver); % Record solver in options.solver; will be used in postprima
% When the invoker is prima, options.solver = solver = '' unless the user defines
% an options.solver in solver_list. Here, '' is an empty char array to signify
% that the solver is yet to decide.

% Check unknown fields according to solver
% solver is '' if it has not been decided yet; in that case, we suppose (for
% simplicity) that all possible fields are known.
known_fields = {'iprint', 'maxfun', 'rhobeg', 'rhoend', 'ftarget', 'classical', 'quiet', 'debug', 'chkfunval', 'solver', 'maxhist', 'output_xhist', 'fortran', 'precision'};
if ~isfield(options, 'classical') || (islogicalscalar(options.classical) && ~options.classical)
    known_fields = [known_fields, 'eta1', 'eta2', 'gamma1', 'gamma2'];
end
if isempty(solver) || any(strcmpi(solver, {'newuoa', 'bobyqa', 'lincoa'}))
    known_fields = [known_fields, 'npt'];
end
if isempty(solver) || any(strcmpi(solver, {'bobyqa', 'lincoa', 'cobyla'}))
    known_fields = [known_fields, 'scale'];
end
if isempty(solver) || strcmpi(solver, 'bobyqa')
    known_fields = [known_fields, 'honour_x0'];
end
if isempty(solver) || any(strcmpi(solver, {'lincoa', 'cobyla'}))
    known_fields = [known_fields, {'ctol', 'cweight', 'maxfilt'}];
end
if isempty(solver) || strcmpi(solver, 'cobyla')
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
        wmsg = sprintf('%s: unknown option %s; it is ignored.', invoker, strjoin(unknown_fields, ', '));
    else
        wmsg = sprintf('%s: unknown options %s; they are ignored.', invoker, strjoin(unknown_fields, ', '));
    end
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Set default npt according to solver
% If solver = '' (empty char array), then invoker must be prima, and a solver
% will be selected later; when the solver is chosen, a valid npt will be defined.
% Note we have to take maxfun into consideration when selecting the solver,
% because npt < maxfun-1 is needed! See function select_solver for details.
if isempty(solver)
    npt = NaN; % The real npt will be (and should be) set when solver is selected
else
    switch lower(solver)
    case {'newuoa', 'bobyqa', 'lincoa'}
        npt = 2*lenx0 + 1;
    case {'uobyqa'}
        npt = (lenx0+1)*(lenx0+2)/2;
    case {'cobyla'}
        npt = lenx0+1;
		% uobyqa and cobyla do not need npt an option, but we need npt to validate/set maxfun
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
if strcmpi(solver, 'bobyqa') && ~options.scale
    rhobeg = max(eps, min(rhobeg, min(ub-lb)/4));
    rhoend = max(eps, min(0.1*rhobeg, rhoend));
end


% Validate the user-specified options; adopt the default values if needed

% Validate options.npt
% There are the following possibilities.
% 1. The user specifies options.npt
% 1.1. The solver is yet to decide (solver = ''): we keep options.npt if it is
% a positive integer; otherwise, raise a warning and set options.npt to NaN;
% 1.2. The user has chosen a valid solver: we keep options.npt if it is
% compatible with the solver; otherwise, raise a warning and set options.npt
% to the default value according to the solver.
% 2. The user does not specify options.npt
% 1.1. The solver is yet to decide (solver = ''): we set options.npt to NaN.
% 1.2. The user has chosen a valid solver: we set options.npt to the default
% value according to the solver.
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
    elseif any(strcmpi(solver, {'newuoa', 'bobyqa', 'lincoa'})) && (~isintegerscalar(options.npt) || isnan(options.npt) || options.npt < lenx0+2 || options.npt > (lenx0+1)*(lenx0+2)/2)
        % newuoa, bobyqa and lincoa requires n+2 <= npt <= (n+1)*)(n+2)/2;
        % uobyqa and cobyla do not use npt.
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
    % When solver = '' (empty char array), the default npt is NaN.
    % For uobyqa and cobyla, we also adopt the 'default npt' defined above,
    % although it will NOT be used by the solver
end
options.npt = double(options.npt);  % All integers will be passed as doubles to the Fortran MEX.
% Although npt and maxfun are integers logically, they have to be passed to the mexified code as
% double variables. In mex, data is passed by pointers, but there are only very limited functions
% that can read an integer value from a pointer or write an integer value to a pointer
% (mxCopyPtrToInteger1, mxCopyInteger1ToPtr, mxCopyPtrToInteger2, mxCopyInteger2ToPtr,
% mxCopyPtrToInteger4, mxCopyInteger4ToPtr; no function for integer*8). This makes it impossible to
% pass integer data properly unless we know the kind of the integer. Therefore, in general, it is
% recommended to pass integers as double variables and then cast them back to integers when needed.
% Indeed, in MATLAB, even if we define npt = 1000, the class of npt is double! To get an integer
% npt, we would have to define npt = int32(1000) or npt = int64(1000)!

% Validate options.maxfun
validated = false;
if isfield(options, 'maxfun')
    if ~isintegerscalar(options.maxfun) || options.maxfun <= 0 || isnan(options.maxfun)
        wid = sprintf('%s:InvalidMaxfun', invoker);
        wmsg = sprintf('%s: invalid maxfun; it should be a positive integer; it is set to %d.', invoker, maxfun);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif options.maxfun > maxint()
        % maxfun would suffer from overflow in the Fortran code
        wid = sprintf('%s:MaxfunTooLarge', funname);
        wmsg = sprintf('%s: maxfun exceeds the upper limit supported; it is set to %d.', funname, maxint());
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
        options.maxfun = maxint();
        validated = true;  % We have set options.maxfun to a valid value in the last line.
    elseif isempty(solver) && options.maxfun <= lenx0+1  % Here, options.maxfun cannot be NaN. No worry about the comparison.
        options.maxfun = lenx0+2; % Here we take lenx0+2 (the smallest possible value for npt)
        validated = true; %!!! % Set validated = true so that options.maxfun will not be set to the default value later
        wid = sprintf('%s:InvalidMaxfun', invoker);
        wmsg = sprintf('%s: invalid maxfun; it should be a positive integer at least n+2; it is set to n+2.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif ~isempty(solver) && options.maxfun <= options.npt  % Here, options.maxfun or options.npt cannot be NaN. No worry about the comparison.
        options.maxfun = options.npt+1; % Here we take npt+1 instead of the default maxfun
        validated = true; %!!! % Set validated = true so that options.maxfun will not be set to the default value later
        wid =  sprintf('%s:InvalidMaxfun', invoker);
        switch lower(solver) % The warning message depends on solver
        case {'newuoa', 'lincoa', 'bobyqa'}
            wmsg = sprintf('%s: invalid maxfun; %s requires maxfun > npt; it is set to npt+1.', invoker, solver);
        case 'uobyqa'
            wmsg = sprintf('%s: invalid maxfun; %s requires maxfun > (n+1)*(n+2)/2; it is set to (n+1)*(n+2)/2+1.', invoker, solver);
        case 'cobyla'
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
options.maxfun = double(options.maxfun);   % All integers will be passed as doubles to the Fortran MEX.
% One can check that options.maxfun >= n+2;

% Validate options.rhobeg
% NOTE: if the problem is to be scaled, then options.rhobeg and options.rhoend
% will be used as the initial and final trust-region radii for the scaled problem.
validated = false;
if isfield(options, 'rhobeg')
    if ~isrealscalar(options.rhobeg) || options.rhobeg <= 0 || isnan(options.rhobeg) || options.rhobeg == inf
        wid = sprintf('%s:InvalidRhobeg', invoker);
        wmsg = sprintf('%s: invalid rhobeg; it should be a positive number; it is set to max(%g, rhoend).', invoker, rhobeg);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif strcmpi(solver, 'bobyqa')  % Validate options.rhobeg for bobyqa
        if options.scale && options.rhobeg > 1  % This case cannot be combined with the next case, as ub and lb are NOT scaled yet in preprima
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
        validated = true; %!!! % Set validated = true so that options.rhobeg will not be set to the default value later
    else
        validated = true;
    end
end
if ~validated % options.rhobeg has not got a valid value yet
    % Take into account `rhoend` if it has got a valid value. We do not do this for `bobyqa`, which
    % requires that rhobeg <= min(xu-xl)/2.
    if isfield(options, 'rhoend') && isrealscalar(options.rhoend) && options.rhoend >= 0 && options.rhoend < inf && ~strcmpi(solver, 'bobyqa')
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
        wmsg = sprintf('%s: invalid rhoend; we should have rhobeg >= rhoend > 0; it is set to min(0.1*rhobeg, %g).', invoker, rhoend);
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
        wmsg = sprintf('%s: invalid ftarget; it should be a real number; it is set to %g.', invoker, ftarget);
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
        wmsg = sprintf('%s: invalid ctol; it should be a nonnegative number; it is set to %g.', invoker, ctol);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated
    options.ctol = ctol;
end
options.ctol = double(options.ctol);  % ctol can be 0

% Validate options.cweight
validated = false;
if isfield(options, 'cweight')
    if ~isrealscalar(options.cweight) || options.cweight < 0 || isnan(options.cweight)
        wid = sprintf('%s:InvalidCweight', invoker);
        wmsg = sprintf('%s: invalid cweight; it should be a nonnegative number; it is set to %g.', invoker, cweight);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated
    options.cweight = cweight;
end
options.cweight = double(options.cweight);  % cweight can be +Inf

% Validate options.classical
validated = false;
if isfield(options, 'classical')
    if ~islogicalscalar(options.classical)
        wid = sprintf('%s:InvalidClassicalFlag', invoker);
        wmsg = sprintf('%s: invalid classical flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(classical));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif options.classical && ~ismember('classical', all_variants())
        wid = sprintf('%s:ClassicalUnavailable', invoker);
        wmsg = sprintf('%s: classical = true but the classical version is unavailable; classical is set to false.', invoker);
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
    wmsg = sprintf('%s: in classical mode, which may CRASH your computer; it is discouraged except for research purposes; set options.classical to false to disable classical mode.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
end

% Validate options.precision
validated = false;
if isfield(options, 'precision')
    if ~ischarstr(options.precision) || ~ismember(lower(options.precision), all_precisions())
        wid = sprintf('%s:InvalidPrecision', invoker);
        wmsg = sprintf('%s: invalid or unavailable precision; it is set to ''%s''.', invoker, precision);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated  % options.precision has not got a validated value yet
    options.precision = precision;
end
options.precision = lower(char(options.precision));

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
        options.fortran = true;
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
        validated = true;
    elseif ~options.fortran && ~strcmp(options.precision, 'double')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % In this version, we support precision ~= 'double' only when calling the Fortran solvers!
        wid = sprintf('%s:FortranContradictPrecision', invoker);
        wmsg = sprintf('%s: fortran = false but precision = %s; fortran is reset to true.', ...
            invoker, options.precision);
        options.fortran = true;
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
        validated = true;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if ~isintegerscalar(options.iprint) || abs(options.iprint) > 3
        wid = sprintf('%s:InvalidIprint', invoker);
        wmsg = sprintf('%s: invalid iprint; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to %d.', invoker, iprint);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
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
    elseif options.iprint ~= 0 && options.classical
        % In classical mode, a nonzero iprint = 0 is not supported.
        wid = sprintf('%s:IprintContradictClassical', invoker);
        wmsg = sprintf('%s: iprint = %d but classical = true; iprint is reset to 0 as other values are not supported.', invoker, options.iprint);
        options.iprint = 0;
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
        validated = true;
    else
        validated = true;
    end
end
if ~validated % options.iprint has not got a valid value yet
    if user_says_quiet
        % The user says "quiet!". Set options.iprint = 0 regardless of the default iprint.
        options.iprint = 0;
    else
        options.iprint = iprint;
    end
end
options.iprint = double(options.iprint);   % All integers will be passed as doubles to the Fortran MEX.

% Validate options.debug
validated = false;
if isfield(options, 'debug')
    if ~islogicalscalar(options.debug)
        wid = sprintf('%s:InvalidDebugflag', invoker);
        wmsg = sprintf('%s: invalid debugging flag; it should be true(1) or false(0); it is set to %s.', invoker, mat2str(debug_flag));
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.debug has not got a valid value yet
    options.debug = debug_flag;
end
options.debug = logical(options.debug);
if options.debug
    wid = sprintf('%s:Debug', invoker);
    wmsg = sprintf('%s: in debug mode; set options.debug to false to disable debug.', invoker);
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
%    if options.quiet
%        options.quiet = false;
%        wid = sprintf('%s:Debug', invoker);
%        wmsg = sprintf('%s: options.quiet is set to false because options.debug = true.', invoker);
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
        wmsg = sprintf('%s: chkfunval = true but debug = false; chkfunval is set to false; set both flags to true to check function values.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif logical(options.chkfunval) && ~strcmp(options.precision, 'double')
        wid = sprintf('%s:InvalidChkfunval', invoker);
        wmsg = sprintf('%s: chkfunval = true but precision = %s; chkfunval is set to false.', invoker, options.precision);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        validated = true;
    end
end
if ~validated % options.chkfunval has not got a valid value yet
    options.chkfunval = logical(chkfunval) && options.debug && strcmp(options.precision, 'double');
end
if options.chkfunval
    wid = sprintf('%s:Chkfunval', invoker);
    if strcmp(solver, 'cobyla')
        wmsg = sprintf('%s: checking whether fx = fun(x) and constr = con(x) at exit, which costs an extra function/constraint evaluation; set options.chkfunval to false to disable the check.', invoker);
    else
        wmsg = sprintf('%s: checking whether fx = fun(x) at exit, which costs an extra function evaluation; set options.chkfunval to false to disable the check.', invoker);
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
        options.maxhist = min(options.maxhist, options.maxfun);
        validated = true;
    end
end
if ~validated  % options.maxhist has not got a valid value
    options.maxhist = options.maxfun;  % options.maxfun has been validated
end
options.maxhist = double(options.maxhist);   % All integers will be passed as doubles to the Fortran MEX.

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
options.output_xhist = logical(options.output_xhist);

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
options.output_nlchist = logical(options.output_nlchist);

% Validate options.maxfilt
validated = false;
if isfield(options, 'maxfilt')
    if ~isintegerscalar(options.maxfilt) || options.maxfilt < 1
        wid = sprintf('%s:InvalidMaxfilt', invoker);
        wmsg = sprintf('%s: invalid maxfilt; it should be a positive integer; it is set to %d.', invoker, maxfilt);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    elseif options.maxfilt < min_maxfilt
        wid = sprintf('%s:InvalidMaxfilt', invoker);
        wmsg = sprintf('%s: maxfilt is too small; it should be an integer at least %d; it is set to %d.', invoker, min_maxfilt, maxfilt);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    else
        options.maxfilt = min(options.maxfilt, options.maxfun);
        validated = true;
    end
end
if ~validated  % options.maxfilt has not got a valid value
    options.maxfilt = maxfilt;  % options.maxfun has been validated
end
options.maxfilt = double(options.maxfilt);  % All integers will be passed as doubles to the Fortran MEX.

% Validate options.eta1
user_eta1_correct = false;  % Does the user provide a correct eta1? Needed when validating eta2.
validated = false;
if isfield(options, 'eta1')
    if ~isrealscalar(options.eta1) || options.eta1 < 0 || options.eta1 >= 1
        wid = sprintf('%s:InvalidEta1', invoker);
        if isfield(options, 'eta2') && isrealscalar(options.eta2) && options.eta2 > 0 && options.eta2 <= 1
        % The user provides a correct eta2; we define eta1 as follows.
            options.eta1 = max(eps, options.eta2/7);
            wmsg = sprintf('%s: invalid eta1; it should be in the interval [0, 1) and not more than eta2; it is set to %g.', invoker, options.eta1);
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
options.eta1 = double(options.eta1);

% Validate options.eta2
validated = false;
if isfield(options, 'eta2')
    if ~isrealscalar(options.eta2) || (isnan(options.eta1) && options.eta2 < 0) || options.eta2 < options.eta1 || options.eta2 > 1
        wid = sprintf('%s:InvalidEta2', invoker);
        if user_eta1_correct
        % The user provides a correct eta1; we define eta2 as follows.
            options.eta2 = (options.eta1 + 2)/3;
            validated = true;
            wmsg = sprintf('%s: invalid eta2; it should be in the interval [0, 1] and not less than eta1; it is set to %g.', invoker, options.eta2);
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
options.eta2 = double(options.eta2);

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
options.gamma1 = double(options.gamma1);

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
options.gamma2 = double(options.gamma2);

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
funname = callstack(1).name; % Name of the current function

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
    % Private/unexpected error
    error(sprintf('%s:InvalidScaling', funname), '%s: UNEXPECTED ERROR: at least one of [-lb; ub] is infinity. Scaling should not be performed.', funname);
end

scaling_factor = (ub - lb)/2;
shift = (lb + ub)/2;

fun = @(x) fun(scaling_factor.*x+shift);
x0 = (x0-shift) ./ scaling_factor;
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
% <==> (lb-shift) ./ scaling_factor < x_after_scaling < (ub-shift) ./ scaling_factor
    lb = (lb-shift) ./ scaling_factor;
end
if ~isempty(ub)
    ub = (ub-shift) ./ scaling_factor;
end
if ~isempty(nonlcon)
    nonlcon = @(x) nonlcon(scaling_factor.*x+shift);
end

% Zaikun 2020-05-25: We do not warn about scaling anymore. Scaling works
% well in several real problems.
%if any(scaling_factor ~= 1)
%    wid = sprintf('%s:ProblemScaled', invoker);
%    wmsg = sprintf('%s: problem scaled according to bound constraints; do this only if the bounds reflect the scaling of variables; if not, set options.scale to false to disable scaling.', invoker);
%    warning(wid, '%s', wmsg);
%    warnings = [warnings, wmsg];
%end

substantially_scaled = false;
%if (max([scaling_factor; 1 ./ scaling_factor]) > substantially_scaled_threshold)
if max([1; scaling_factor])/min([1; scaling_factor]) > substantially_scaled_threshold
    substantially_scaled = true;
end

if min(scaling_factor) < eps
    % Private/unexpected error
    error(sprintf('%s:InvalidScaling', funname), '%s: UNEXPECTED ERROR: invalid scaling factor returned when called by %s.', funname, invoker);
end
return

%%%%%%%%%%%%%%%%%%%%%%%% Function for selecting solver %%%%%%%%%%%%%%%%%%%%
function [options, warnings] = select_solver(invoker, options, probinfo, warnings)

invoker_list = {'prima'};
% Only prima needs select_solver. We may have other invokers in the future!
solver_list = all_solvers();
% We may add other solvers in the future!
% Note that prima is not a possible solver here!
callstack = dbstack;
funname = callstack(1).name; % Name of the current function

if ~ismember(invoker, invoker_list)
    % Private/unexpected error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s serves only %s.', funname, funname, strjoin(invoker_list, ', '));
end
% After pre_options, options.solver is either a member of solver_list
% or '' (i.e., an empty char array), the second signifying the solver
% is yet to decide.
% 1. If options.solver is in solver_list, we check whether it can solve the
% problem. If yes, we set solver = options.solver; otherwise, we warn about
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
            solver = 'uobyqa';
        elseif (options.maxfun <= n+2) % After preprima, options.maxfun >= n+2 is ensured. Thus options.maxfun <= n+2 indeed means options.maxfun = n+2
            solver = 'cobyla';
        else
            solver = 'newuoa';  % options.npt will be set later
            % Interestingly, we note in our test that lincoa outperformed
            % newuoa on unconstrained CUTEst problems when the dimension
            % was not large (i.e., <= 50) or the precision requirement
            % was not high (i.e., >= 1e-5). Therefore, it is worthwhile to
            % try lincoa when an unconstrained problem is given.
            % Nevertheless, for the moment, we set the default solver
            % for unconstrained problems to be newuoa.
        end
    case 'bound-constrained'
        if (options.maxfun <= n+2)
            solver = 'cobyla';
        else
            solver = 'bobyqa';  % options.npt will be set later
        end
    case 'linearly-constrained'
        if (options.maxfun <= n+2)
            solver = 'cobyla';
        else
            solver = 'lincoa';  % options.npt will be set later
        end
    case 'nonlinearly-constrained'
        solver = 'cobyla';
    otherwise
        % Private/unexpected error
        error(sprintf('%s:InvalidProbType', funname), '%s: UNEXPECTED ERROR: unknown problem type ''%s'' received.', funname, ptype);
    end
end

% Revise options.npt according to the selected solver
% Note that pre_options has set options.npt to either a positive integer or NaN.
if ismember(solver, {'newuoa', 'bobyqa', 'lincoa'}) && (isnan(options.npt) || options.npt < n+2 || options.npt > min((n+1)*(n+2)/2, options.maxfun-1))
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
if strcmp(solver, 'bobyqa') && options.rhobeg > min(probinfo.refined_data.ub-probinfo.refined_data.lb)/2
    options.rhobeg = max(eps, min(probinfo.refined_data.ub-probinfo.refined_data.lb)/4);
    options.rhoend = max(eps, min(0.1*options.rhobeg, options.rhoend));
    if ismember('rhobeg', probinfo.user_options_fields) || ismember('rhoend', probinfo.user_options_fields)
        wid = sprintf('%s:InvalidRhobeg', invoker);
        wmsg = sprintf('%s: rhobeg is set to %g and rhoend to %g according to the selected solver bobyqa, which requires rhoend <= rhobeg <= min(ub-lb)/2.', invoker, options.rhobeg, options.rhoend);
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

solver_list = [all_solvers(), 'prima'];
% Note: prima is also a possible solver when prob_solv_match is called in
% preprima to check whether the invoker can handle the problem.

if ~ismember(solver, solver_list)
    % Private/unexpected error
    error(sprintf('%s:InvalidSolver', funname), ...
    '%s: UNEXPECTED ERROR: unknown solver ''%s'' received.', funname, solver);
end

switch ptype
case 'unconstrained'
    match = true;
    % Essentially do nothing. DO NOT remove this case. Otherwise, the
    % case would be included in 'otherwise', which is not correct.
case 'bound-constrained'
    match = ismember(solver, all_solvers('bound_constrained_solvers'));
case 'linearly-constrained'
    match = ismember(solver, all_solvers('linearly_constrained_solvers'));
case 'nonlinearly-constrained'
    match = ismember(solver, all_solvers('nonlinearly_constrained_solvers'));
otherwise
    % Private/unexpected error
    error(sprintf('%s:InvalidProbType', funname), '%s: UNEXPECTED ERROR: unknown problem type ''%s'' received.', funname, ptype);
end

match = match || strcmp(solver, 'prima');  % prima matches all types of problems.

return

%%%%%%%%%%%%%%% Function for calculating constraint violation %%%%%%%%%%
function [constrviolation, nlcineq, nlceq] = get_constrv(x, Aineq, bineq, Aeq, beq, lb, ub, nonlcon)
%GET_CONSTRV calculates the absolute constraint violation at x.
nlcineq = [];
nlceq = [];
if ~isempty(nonlcon)
    [nlcineq, nlceq] = nonlcon(x);
end
constrviolation = get_cstrv(x, Aineq, bineq, Aeq, beq, lb, ub, nlcineq, nlceq);
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

solver_list = {'bobyqa'}; % Only BOBYQA needs pre_rhobeg_x0. May have others in the future.

if ~ismember(lower(options.solver), solver_list)
    % Private/unexpected error
    error(sprintf('%s:InvalidSolver', funname), '%s: UNEXPECTED ERROR: %s serves only %s.', funname, funname, strjoin(solver_list, ', '));
end

if isfield(options, 'honour_x0') && options.honour_x0  % In this case, we respect the user-defined x0 and revise rhobeg
    rhobeg_old = options.rhobeg;
    lbx = (lb > -inf & x0 - lb <= eps*max(abs(lb), 1));  % x0 essentially equals lb
    ubx = (ub < inf & x0 - ub >= -eps*max(abs(ub), 1));  % x0 essentially equals ub
    x0(lbx) = lb(lbx);
    x0(ubx) = ub(ubx);
    options.rhobeg = max(eps, min([options.rhobeg; x0(~lbx) - lb(~lbx); ub(~ubx) - x0(~ubx)]));
    if rhobeg_old - options.rhobeg > eps*max(1, rhobeg_old)
        options.rhoend = max(eps, min(0.1*options.rhobeg, options.rhoend));  % We do not revise rhoend unless rhobeg is revised
        if ismember('rhobeg', user_options_fields) || ismember('rhoend', user_options_fields)
            wid = sprintf('%s:ReviseRhobeg', invoker);
            wmsg = sprintf('%s: rhobeg is revised to %g and rhoend to %g so that the distance between x0 and the inactive bounds is at least rhobeg.', invoker, options.rhobeg, options.rhoend);
            warning(wid, '%s', wmsg);
            warnings = [warnings, wmsg];
        end
    else
        options.rhoend = min(options.rhoend, options.rhobeg);  % This may update rhoend slightly
    end
else
    % N.B.: The following code is valid only if lb <= x0 <= ub and rhobeg <= min(ub-lb)/2, which
    % hold after `pre_options` and `project` are invoked.
    x0_old = x0;
    lbx = (x0 <= lb + 0.5*options.rhobeg);
    lbx_plus = (x0 > lb + 0.5*options.rhobeg) & (x0 < lb + options.rhobeg);
    ubx_minus = (x0 < ub - 0.5*options.rhobeg) & (x0 > ub - options.rhobeg);
    ubx = (x0 >= ub - 0.5*options.rhobeg);
    x0(lbx) = lb(lbx);
    x0(lbx_plus) = lb(lbx_plus) + options.rhobeg;
    x0(ubx_minus) = ub(ubx_minus) - options.rhobeg;
    x0(ubx) = ub(ubx);
    if any(abs(x0_old-x0) > 0)
        wid = sprintf('%s:ReviseX0', invoker);
        wmsg = sprintf('%s: x0 is revised so that the distance between x0 and the inactive bounds is at least rhobeg; set options.honour_x0 to true if you prefer to keep x0.', invoker);
        warning(wid, '%s', wmsg);
        warnings = [warnings, wmsg];
    end
end
return


%% Function for checking whether the Fortran MEX function is available and revising options.fortran %
function [options, warnings] = is_fortran_available(invoker, options, warnings)
% IS_FORTRAN_AVAILABLE checks whether the Fortran MEX function is available. If no, raise a warning
% and set options.fortran to false. It does nothing if options.fortran is false or does not exist.
if ~isfield(options, 'fortran') || ~options.fortran
    return
end
mfiledir = fileparts(mfilename('fullpath'));  % The directory where this .m file resides.
mexdir = mfiledir;  % The directory where the MEX functions reside.
if options.classical
    variant = 'classical';
else
    variant = 'modern';
end
mexname = get_mexname(options.solver, options.precision, options.debug, variant, mexdir);
if exist(fullfile(mexdir, [mexname, '.', mexext]), 'file') ~= 3
    wid = sprintf('%s:FortranNotAvailable', invoker);
    wmsg = sprintf('%s: fortran = true but the Fortran MEX function is not available for the %s variant of %s with precision %s; fortran is reset to false.', invoker, variant, options.solver, options.precision);
    options.fortran = false;
    warning(wid, '%s', wmsg);
    warnings = [warnings, wmsg];
    options.fortran = false;
end
return
