function [x, fx, exitflag, output] = postpdfo(probinfo, output)
%POSTPDFO postprocesses the output by pdfo or its solvers and creates the
%   output variables.
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
% 1. All errors in this function are unexpcted errors, which means they
% should not occur unless there is a bug in the code.
% 2. Some unexpcted errors are external.
%
% TODO: None
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% postpdfo starts

% Obligatory fields in output
% If a new solver is included, it should include at least the following
% fields in output. For unconstrained problems, put constrviolation = 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obligatory_output_fields = {'x', 'fx', 'exitflag', 'funcCount', 'constrviolation'};%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obligatory fields in probinfo and options
obligatory_probinfo_fields = {'raw_data', 'refined_data', 'fixedx', 'fixedx_value', ...
    'nofreex', 'infeasible_bound', 'infeasible_lineq', 'infeasible_leq', ...
    'trivial_lineq', 'trivial_leq', 'infeasible', 'scaled', 'scaling_factor', ...
    'shift', 'reduced', 'raw_type', 'raw_dim', 'refined_type', 'refined_dim', ...
    'feasibility_problem', 'user_options_fields', 'options', 'warnings'};
obligatory_options_fields = {'classical', 'debug', 'chkfunval', 'quiet'};

% All possible solvers
unconstrained_solver_list = {'uobyqa', 'newuoa'};
constrained_solver_list = {'bobyqa', 'lincoa', 'cobyla'};
nonlinearly_constrained_solver_list = {'cobyla'};
solver_list = [unconstrained_solver_list, constrained_solver_list];

% Solvers from PDFO; there may be external solvers added later.
internal_solver_list = {'uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla'};

% Who is calling this function? Is it a correct invoker?
invoker_list = ['pdfo', solver_list];
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if (length(callstack) == 1) || ~ismember(callstack(2).name, invoker_list)
    % Private/unexpcted error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s should only be called by %s.', funname, funname, mystrjoin(invoker_list, ', '));
else
    invoker = callstack(2).name; % Name of the function who calls this function
end

% With extreme barrier (implemented when options.classical=false), all
% the function values that are NaN or larger than hugefun are replaced
% by hugefun; all the constraint values that are NaN or larger than
% hugecon are replaced by hugecon. hugefun and hugecon are defined in
% pdfoconst.F, and can be obtained by gethuge.
hugefun = gethuge('fun');
hugecon = gethuge('con');

% Verify the input before starting the real business
% Verify probinfo
if (length(callstack) >= 3) && strcmp(callstack(3).name, 'pdfo')
% In this case, prepdfo sets probinfo to empty.
    if ~isempty(probinfo)
        % Public/unexpected error
        error(sprintf('%s:InvalidProbinfo', invoker),...
            '%s: UNEXPECTED ERROR: probinfo should be empty because %s is a solver called by pdfo.', invoker, invoker);
    end
else
    if ~isa(probinfo, 'struct')
        % Public/unexpected error
        error(sprintf('%s:InvalidProbinfo', invoker),...
            '%s: UNEXPECTED ERROR: probinfo should be a structure.', invoker);
    end
    missing_fields = setdiff(obligatory_probinfo_fields, fieldnames(probinfo));
    if ~isempty(missing_fields)
        % Public/unexpected error
        error(sprintf('%s:InvalidProbinfo', invoker),...
            '%s: UNEXPECTED ERROR: probinfo misses the %s field(s).', invoker, mystrjoin(missing_fields, ', '));
    end

    % Read and verify options
    options = probinfo.options;
    if ~isa(options, 'struct')
        % Public/unexpected error
        error(sprintf('%s:InvalidOptions', invoker), ...
            '%s: UNEXPECTED ERROR: options should be a structure.', invoker);
    end
    missing_fields = setdiff(obligatory_options_fields, fieldnames(options));
    if ~isempty(missing_fields)
        % Public/unexpected error
        error(sprintf('%s:InvalidOptions', invoker),...
            '%s: UNEXPECTED ERROR: options misses the %s field(s).', invoker, mystrjoin(missing_fields, ', '));
    end
end

% The solver that did the computation (needed for verifying output below)
if strcmp(invoker, 'pdfo')
    % In this case, the invoker is pdfo rather than a solver called by pdfo.
    % Thus probinfo is nonempty, and options has been read and verified as above.
    solver = options.solver;
else
    solver = invoker;
end
if isempty(solver) || (~isa(solver, 'char') && ~isa(solver, 'string')) || ~ismember(solver, solver_list)
    % Public/unexpected error
    error(sprintf('%s:InvalidSolver', invoker), '%s: UNEXPECTED ERROR: invalid solver passed to %s.', invoker, funname);
end

% Verify output
if ~isa(output, 'struct')
    % Public/unexpcted error
    error(sprintf('%s:InvalidOutput', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an output that is not a structure', invoker, solver);
end
if ismember(solver, internal_solver_list)
    % For internal solvers, output should contain fhist, chist, and warnings
    obligatory_output_fields = [obligatory_output_fields, 'fhist', 'chist', 'warnings'];
end
if strcmp(solver, 'lincoa')
    % For lincoa, output should contain constr_modified
    obligatory_output_fields = [obligatory_output_fields, 'constr_modified'];
end
if ismember(solver, nonlinearly_constrained_solver_list) && ismember(solver, internal_solver_list)
    % For nonlinearly constrained internal solvers, output should contain nlinceq and nlceq
    obligatory_output_fields = [obligatory_output_fields, 'nlcineq', 'nlceq'];
end
missing_fields = setdiff(obligatory_output_fields, fieldnames(output));
if ~isempty(missing_fields)
    % Public/unexpected error
    error(sprintf('%s:InvalidOutput', invoker),...
        '%s: UNEXPECTED ERROR: %s returns an output that misses the %s field(s).', invoker, solver, mystrjoin(missing_fields, ', '));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the invoker is a solver called by pdfo, then let pdfo do the postprecessing
% Put this after verifying output, because we will use the information in it.
if (length(callstack) >= 3) && strcmp(callstack(3).name, 'pdfo')
    x = output.x;
    fx = output.fx;
    exitflag = output.exitflag;
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Record solver name in output (should be done after verifying that
% output is a structure.
output.algorithm = solver;

% Read information in output
x = output.x;
output = rmfield(output, 'x'); % output does not include x at return
fx = output.fx;
output = rmfield(output, 'fx'); % output does not include fx at return
exitflag = output.exitflag;
output = rmfield(output, 'exitflag'); % output does not include exitflag at return
nf = output.funcCount;
constrviolation = output.constrviolation;
if strcmp(solver, 'lincoa')
    constr_modified = output.constr_modified;
    output = rmfield(output, 'constr_modified');
end
if ~isfield(output, 'warnings') || isempty(output.warnings)
    output.warnings = {};
end

% Verify x
if ~isnumeric(x) || ~isreal(x) || ~isvector(x) || size(x,2)~=1
    % Public/unexpcted error
    error(sprintf('%s:InvalidX', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an x that is not a real column or scalar.', invoker, solver);
end

% Verify fx
if ~isnumeric(fx) || ~isreal(fx) || ~isscalar(fx)
    % Public/unexpcted error
    error(sprintf('%s:InvalidFx', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an fx that is not a real number.', invoker, solver);
end

% Verify exitflag
if ~isnumeric(exitflag) || ~isscalar(exitflag) || ~isreal(exitflag) || rem(exitflag, 1)~=0
    % Public/unexpcted error
    error(sprintf('%s:InvalidExitFlag', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an exitflag that is not an integer', invoker, solver);
end

% Verify nf
if ~isnumeric(nf) || ~isscalar(nf) || ~isreal(nf) || rem(nf, 1)~=0 || nf < 0
    % Public/unexpcted error
    error(sprintf('%s:InvalidNF', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an nf that is not a nonnegative integer.', invoker, solver);
end
if nf <= 0
    % If prepdfo works properly, then nf<=0 should never happen.
    % Public/unexpcted error
    error(sprintf('%s:InvalidNF', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns nf=0 unexpectedly with exitflag %d.', invoker, solver, exitflag);
end

% Read and verify fhist
if isfield(output, 'fhist')
    fhist = output.fhist;
else % External solvers may not return fhist
    fhist = fx + zeros(1, nf);
end
if ~isempty(fhist) && (~isnumeric(fhist) || ~isreal(fhist) || ~isvector(fhist) || (nf ~= length(fhist)))
    % Public/unexpected error
    error(sprintf('%s:InvalidFhist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an fhist that is not a real vector of length nf.', invoker, solver);
end
if ~options.classical && ~probinfo.infeasible && ~probinfo.nofreex
    if any(fhist > hugefun) || any(isnan(fhist))
        % Public/unexpected error
        error(sprintf('%s:InvalidFhist', invoker), ...
             '%s: UNEXPECTED ERROR: %s returns an fhist with NaN or values larger than hugefun=%1.2e; this is impossible with extreme barrier.', invoker, solver, hugefun);
    elseif ~isempty(fhist) && max(fhist) == hugefun
        wid = sprintf('%s:ExtremeBarrier', invoker);
        wmessage = sprintf('%s: extreme barrier is invoked; function values that are NaN or larger than hugefun=%1.2e are replaced by hugefun.', invoker, hugefun);
        warning(wid, '%s', wmessage);
        output.warnings = [output.warnings, wmessage];
    end
end

% If the problem is a feasibility problem, set fx to [], and remove fhist from output;
if probinfo.feasibility_problem
    fx = [];
    output = rmfield(output, 'fhist');
    if ~strcmp(probinfo.refined_type, 'nonlinearly-constrained')
        % No function evaluation involved when solving a linear feasibility problem.
        % By "function evaluation", we mean the evaluation of the objective function
        % and nonlinear constraint functions, which do not exist in this case.
        % For nonlinear feasibility problems, funcCount is positive.
        output.funcCount = 0;
    end
end

% Verify constrviolation
if ~isnumeric(constrviolation) || ~isscalar(constrviolation) || ~isreal(constrviolation)
    % Public/unexpected error
    error(sprintf('%s:InvalidConstrViolation', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns a constrviolation that is not a real number.', invoker, solver)
end

% Read and verify chist
if isfield(output, 'chist')
    output_has_chist = true;
    chist = output.chist;
else % External solvers may not return chist
    output_has_chist = false;
    chist = constrviolation + zeros(1, nf);
end
if ~(isempty(chist) && ismember(solver, unconstrained_solver_list)) && (~isnumeric(chist) || ~isreal(chist) || ~isvector(chist) || (nf ~= length(chist)))
    % Public/unexpected error
    error(sprintf('%s:InvalidChist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns a chist that is not a real vector of length nf.', invoker, solver);
end
if ~options.classical && ~probinfo.infeasible && ~probinfo.nofreex
    if strcmp(solver, 'cobyla') && (any(chist > hugecon) || any(isnan(chist)))
        % Public/unexpected error
        error(sprintf('%s:InvalidChist', invoker), ...
             '%s: UNEXPECTED ERROR: %s returns a chist with NaN or values larger than hugecon=%1.2e; this is impossible with extreme barrier.', invoker, solver, hugecon);
    elseif ~isempty(chist) && max(chist) == hugecon
        wid = sprintf('%s:ExtremeBarrier', invoker);
        wmessage = sprintf('%s: extreme barrier is invoked; constraint values that are NaN or larger than hugecon=%1.2e are replaced by hugecon.', invoker, hugecon);
        warning(wid, '%s', wmessage);
        output.warnings = [output.warnings, wmessage];
    end
end

% Read and verify nlcineq and nlceq
if isfield(output, 'nlcineq')
    output_has_nlcineq = true;
    nlcineq = output.nlcineq;
else
    output_has_nlcineq = false;
    nlcineq = [];
end
if isfield(output, 'nlceq')
    output_has_nlceq = true;
    nlceq = output.nlceq;
else
    output_has_nlceq = false;
    nlceq = [];
end
if ~strcmp(probinfo.refined_type, 'nonlinearly-constrained') && (isfield(output, 'nlcineq') && ~isempty(output.nlcineq) || isfield(output, 'nlceq') && ~isempty(output.nlceq))
    % Public/unexpected error
    error(sprintf('%s:InvalidNonlinearConstraint', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns values of nonlinear constraints for a problem without such constraints.', invoker, solver);
end
if (isfield(output, 'nlcineq') && ~isfield(output, 'nlceq')) || (~isfield(output, 'nlcineq') && isfield(output, 'nlceq'))
    % Public/unexpected error
    error(sprintf('%s:InvalidNonlinearConstraint', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns only one of nlcineq and nlceq; it should return both of them or neither of them.', invoker, solver);
end

% After verification, extract and process the data.

% The problem was (possibly) scaled. Scale it back.
% The scaling affects constrviolation when there are bound constraint.
% Hence constrviolation has to be recalculated so that it equals the
% constraint violation of the returned x with respect to the original problem.
% Ideally, chist should also be recalculated. However, it is impossible
% because we do not save the history of x. Therefore, when
% probinfo.scaled=true, chist is not the history of constraint violation
% of the original problem but the scaled one. It it not consistent with
% constrviolation. Without saving of history of x, we cannot do better.
%
% Before recalculating constrviolation, save the one returned by the
% solver, because it will be used in debug mode when checking whether fx
% is consistent with fhist and chist. See the definition of fhistf for
% details.
constrv_returned = constrviolation;
if probinfo.scaled
    % First calculate the residuals of the linear constraints. This must
    % be calculated before x is scaled back. Otherwise, we would have to
    % scale also the linear constraints back to get the correct residuals.
    % Note that we cannot use probinfo.raw_data, which is available only
    % in debug mode.
    Aineq = probinfo.refined_data.Aineq;
    bineq = probinfo.refined_data.bineq;
    Aeq = probinfo.refined_data.Aeq;
    beq = probinfo.refined_data.beq;
    rineq = [];
    req = [];
    if ~isempty(Aineq)
        rineq = Aineq*x-bineq;
    end
    if ~isempty(Aeq)
        req = Aeq*x-beq;
    end

    % Scale x back
    x = probinfo.scaling_factor.*x + probinfo.shift;

    % Scale bounds back
    lb = probinfo.scaling_factor.*probinfo.refined_data.lb + probinfo.shift;
    ub = probinfo.scaling_factor.*probinfo.refined_data.ub + probinfo.shift;
    if isempty(lb)
        lb = -inf(size(x));
    end
    if isempty(ub)
        ub = inf(size(x));
    end

    % We only need to calculate constrviolation for lincoa and cobyla,
    % because uobyqa and newuoa do not handle constrained problems,
    % while bobyqa is a feasible method and should return constrviolation=0
    % regardless of the scaling unless something goes wrong.
    if strcmp(solver, 'lincoa')
        constrviolation = max([0; rineq; abs(req); lb-x; x-ub], [], 'includenan');
        % max(X, [], 'includenan') returns NaN if X contains NaN, and
        % maximum of X otherwise
    else
        constrviolation = max([0; rineq; abs(req); lb-x; x-ub; nlcineq; abs(nlceq)], [], 'includenan');
        % max(X, [], 'includenan') returns NaN if X contains NaN, and
        % maximum of X otherwise
    end
end

% The problem was (possibly) reduced. Get the full x.
if probinfo.reduced
    freex_value = x;
    x = NaN(length(x)+length(probinfo.fixedx_value), 1);
    x(probinfo.fixedx) = probinfo.fixedx_value;
    x(~probinfo.fixedx) = freex_value;
end

% Set output.constrviolation to the revised constraint violation
output.constrviolation = constrviolation;

% Revise output.constrviolation and output.chist according to problem type
if strcmp(probinfo.refined_type, 'unconstrained') && (constrviolation > 0 || max([0, chist]) > 0)
    % Public/unexpected error
    error(sprintf('%s:InvalidConstrViolation', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns positive constrviolations for an unconstrained problem.', invoker, solver);
end
if strcmp(probinfo.raw_type, 'unconstrained')
    % Do not include constrviolation or chist in output for unconstrained problems
    if isfield(output, 'constrviolation')
        output = rmfield(output, 'constrviolation');
    end
    if isfield(output, 'chist')
        output = rmfield(output, 'chist');
    end
elseif strcmp(probinfo.refined_type, 'unconstrained') && ~strcmp(probinfo.raw_type, 'unconstrained')
    output.constrviolation = 0.0;
    if output_has_chist
        output.chist = zeros(1, nf);
    end
end

% Revise output.nlcineq and output.nlceq according to problem type
if ~strcmp(probinfo.raw_type, 'nonlinearly-constrained')
    if isfield(output, 'nlcineq')
        output = rmfield(output, 'nlcineq');
    end
    if isfield(output, 'nlceq')
        output = rmfield(output, 'nlceq');
    end
end

% Record the retrun message in output.message according to exitflag
switch exitflag % If prepdfo works properly, then 5, 6, 10, 11, 12 should never happen
case 0
    output.message = sprintf('Return from %s because the trust region radius reaches its lower bound.', solver);
case 1
    output.message = sprintf('Return from %s because the target function value is achieved.', solver);
case 2
    output.message = sprintf('Return from %s because a trust region step has failed to reduce the quadratic model.', solver);
case 3
    output.message = sprintf('Return from %s because the objective function has been evaluated maxfun times.', solver);
case 4
    output.message = sprintf('Return from %s because of much cancellation in a denominator.', solver);
%case 5
%    output.message = sprintf('Return from %s because npt is not in the required interval.', solver);
%case 6
%    output.message = sprintf('Return from %s because one of the differences xu(i) - xl(i) is less than 2*rhobeg.', solver);
case 7
    output.message = sprintf('Return from %s because rounding errors are becoming damaging.', solver);
case 8
    output.message = sprintf('Return from %s because rounding errors prevent reasonable changes to x.', solver);
case 9
    output.message = sprintf('Return from %s because the denominator of the updating formula is zero.', solver);
%case 10
%    output.message = sprintf('Return from %s because n should not be less than 2.', solver);
%case 11
%    output.message = sprintf('Return from %s because maxfun is less than npt+1.', solver);
%case 12
%    output.message = sprintf('Return from %s because the gradient of a constraint is zero.', solver);
case 13
    output.message = sprintf('Return from %s because all the variables are fixed by the bounds.', invoker);
case 14
    output.message = sprintf('%s receives a linear feasibility problem and finds a feasible point.', invoker);
case 15
    output.message = sprintf('%s receives a linear feasibility problem but does not find a feasible point.', invoker);
case -1
    output.message = sprintf('Return from %s because NaN occurs in x.', solver);
case -2
    if strcmp(solver, 'cobyla')
        output.message = sprintf('Return from %s because the objective function returns an NaN or nearly infinite value, or the constraints return a NaN.', solver);
    else
        output.message = sprintf('Return from %s because the objective function returns an NaN or nearly infinite value.', solver);
    end
case -3
    output.message = sprintf('Return from %s because NaN occurs in the models.', solver);
case -4
    % Record indices of infeasible constraints
    if any(probinfo.infeasible_lineq)
        output.InfeasibleLinearIneq = find(probinfo.infeasible_lineq)';
        % 'find' changes an vector of true/false to a vector containing the indixes of the true values
    end
    if any(probinfo.infeasible_leq)
        output.InfeasibleLinearEq = find(probinfo.infeasible_leq)';
    end
    if any(probinfo.infeasible_bound)
        output.InfeasibleBound = find(probinfo.infeasible_bound)';
    end
    output.message = sprintf('Return from %s because the constraints are infeasible.', invoker);
otherwise
    % Public/unexpected error
    error(sprintf('%s:InvalidExitflag', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an invalid exitflag %d.', invoker, solver, exitflag);
end

% Record indices of trivial constraints
if any(probinfo.trivial_lineq)
    output.TrivialLinearIneq = find(probinfo.trivial_lineq)';
end
if any(probinfo.trivial_leq)
    output.TrivialLinearEq = find(probinfo.trivial_leq)';
end

% Record warnings in output.warnings
if isfield(output, 'warnings')
    warnings = output.warnings;
    output = rmfield(output, 'warnings');
    % warnings is removed from output and rejoined later, so that it will be the last element of output
else
    warnings = {};
end
if isfield(probinfo, 'warnings')
    warnings = [probinfo.warnings, warnings];
end
if ~isempty(warnings)
    output.warnings = warnings;
%    if ~options.quiet
%        fprintf('The following warnings were rasied by %s:\n', invoker);
%        fprintf('%s\n', warnings{:});
%    end
end

% Print the return message if options.quiet=false
if ~options.quiet
    fprintf('\n%s\n\n', output.message)
end

% Recover the default warning behavior of displaying stack trace, which was disabled by pdfo or its solvers
warning('on', 'backtrace');

% At this point, we have completed defining the outputs (i.e., x, fx,
% exitflag, and output). They will NOT (should NOT) be revised any more.
% The remaining code is reachable only in debug mode.


% More careful checks about fx, constrviolation, fhist, and chist.
% We do this only if the code is in debug mode but not in classical
% mode. The classical mode cannot pass these checks.
if options.debug && ~options.classical
    % Check whether fx is 'optimal'
    fhistf = fhist;
    if ismember(solver, constrained_solver_list)
        fhistf = fhistf(chist <= max(constrv_returned, 0));
    end
    minf = min([fhistf,fx]);
%% Zaikun 2021-05-26: The following test is disabled for lincoa for the moment. lincoa may not pass it.
%%    if (fx ~= minf) && ~(isnan(fx) && isnan(minf)) && ~(strcmp(solver, 'lincoa') && constr_modified)
    if (fx ~= minf) && ~(isnan(fx) && isnan(minf)) && ~strcmp(solver, 'lincoa')
        % Public/unexpected error
        error(sprintf('%s:InvalidFhist', invoker), ...
             '%s: UNEXPECTED ERROR: %s returns an fhist that does not match nf or fx.', invoker, solver);
    end

    % Check whether constrviolation is correct
    cobyla_prec = 1e-6;
    lincoa_prec = 1e-12;
    % COBYLA cannot ensure fx=fun(x) or conval=con(x) due to rounding
    % errors. Instead of checking the equality, we check whether the
    % relative error is within cobyla_prec.
    % There can also be a difference between constrviolation and conv due
    % to rounding errors, especially if the problem is scaled.
    constrviolation = 0;
    if isfield(output, 'constrviolation')
        constrviolation = output.constrviolation;
    end
    if strcmp(solver, 'bobyqa') && (max([chist, constrviolation]) > 0) && ~probinfo.infeasible
        % Public/unexpected error
        error(sprintf('%s:InvalidChist', invoker), ...
             '%s: UNEXPECTED ERROR: %s is a feasible solver yet it returns positive constrviolations.', invoker, solver);
    end
    if (strcmp(solver, 'lincoa') && ~constr_modified) || strcmp(solver, 'cobyla')
        Aineq = probinfo.raw_data.Aineq;
        bineq = probinfo.raw_data.bineq;
        Aeq = probinfo.raw_data.Aeq;
        beq = probinfo.raw_data.beq;
        lb = probinfo.raw_data.lb(:);
        ub = probinfo.raw_data.ub(:);
        lb(isnan(lb)) = -inf; % Replace the NaN in lb with -inf
        ub(isnan(ub)) = inf; % Replace the NaN in ub with inf
        bineq(isnan(bineq)) = inf; % Replace the NaN in bineq with inf
        if ~isempty(Aeq)
            nan_eq = isnan(sum(abs(Aeq), 2)) & isnan(beq); % NaN equality constraints
            Aeq = Aeq(~nan_eq, :); % Remove NaN equality constraints
            beq = beq(~nan_eq);
        end
        if isempty(lb)
            lb = -inf(size(x));
        end
        if isempty(ub)
            ub = inf(size(x));
        end
        rineq = [];
        req = [];
        if ~isempty(Aineq)
            rineq = Aineq*x-bineq;
        end
        if ~isempty(Aeq)
            req = Aeq*x-beq;
        end
        if strcmp(solver, 'lincoa')
            conv = max([0; rineq; abs(req); lb-x; x-ub], [], 'includenan');
            % max(X, [], 'includenan') returns NaN if X contains NaN, and
            % maximum of X otherwise
        else
            conv = max([0; rineq; abs(req); lb-x; x-ub; nlcineq; abs(nlceq)], [], 'includenan');
            % max(X, [], 'includenan') returns NaN if X contains NaN, and
            % maximum of X otherwise
        end

        if ~(isnan(conv) && isnan(constrviolation)) && ~(conv == inf && constrviolation == inf) && ~(abs(constrviolation-conv) <= lincoa_prec*max(1,abs(conv)) && strcmp(solver, 'lincoa')) && ~(abs(constrviolation-conv) <= cobyla_prec*max(1,abs(conv)) && strcmp(solver, 'cobyla'))
            % Public/unexpected error
            error(sprintf('%s:InvalidChist', invoker), ...
              '%s: UNEXPECTED ERROR: %s returns a constrviolation that does not match x.', invoker, solver);
        end
        if isnan(fx)
            cf = chist(isnan(fhist));
        else
            cf = chist(fhist == fx);
        end
        if ~any(cf == constrv_returned) && ~(isnan(constrv_returned) && ~any(~isnan(cf)))
            % Public/unexpected error
            error(sprintf('%s:InvalidFhist', invoker), ...
              '%s: UNEXPECTED ERROR: %s returns a constrviolation that does not match chist.', invoker, solver);
        end
    end

    if options.chkfunval % Check the values of fun(x) and con(x)
        % Check whether fx = fun(x)
        % Recall that probinfo.raw_dat.objective was raw data.
        % When the code arrives here, options.raw_data.objective passed the
        % validation but not preprocessed. It can be empty, a function handle,
        % or a function name. If it is empty, then the objective function used
        % in computation was 0; if it is a function name, then calling it by
        % writing 'objective(x)' will cause an error.
        objective = probinfo.raw_data.objective;
        if isempty(objective)
            funx = 0;
        else
            funx = feval(objective, x);
        end
        if (funx ~= funx) || (funx > hugefun)
            funx = hugefun;
            % Due to extreme barrier (implemented when options.classical=false),
            % all the function values that are NaN or larger than
            % hugefun are replaced by hugefun.
        end
        %if (funx ~= fx) && ~(isnan(fx) && isnan(funx))
        % it seems that COBYLA can return fx~=fun(x) due to rounding
        % errors. Therefore, we cannot use "fx~=funx" to check COBYLA
        if ~(isnan(fx) && isnan(funx)) && ~((fx==funx) || (abs(funx-fx) <= cobyla_prec*max(1, abs(fx)) && strcmp(solver, 'cobyla')))
            % Public/unexpected error
            error(sprintf('%s:InvalidFx', invoker), ...
                '%s: UNEXPECTED ERROR: %s returns an fx that does not match x.', invoker, solver);
        end

        % Check whether [output.nlcineq,  output.nlceq] = nonlcon(x)
        if output_has_nlcineq && output_has_nlceq
            nlcineqx = [];
            nlceqx = [];
            nonlcon = probinfo.raw_data.nonlcon;
            if ~isempty(nonlcon)
                [nlcineqx, nlceqx] = feval(nonlcon, x);
                % Due to extreme barrier (implemented when options.classical=false),
                % all the constraint values that are NaN or above hugecon
                % are replaced by hugecon, and all those below -hugecon are
                % replaced by -hugecon.
                nlcineqx(nlcineqx ~= nlcineqx | nlcineqx > hugecon) = hugecon;
                nlcineqx(nlcineqx < -hugecon) = -hugecon; % This one was done not as an extreme barrier but for numerical stability.
                nlceqx(nlceqx ~= nlceqx | nlceqx > hugecon) = hugecon;
                nlceqx(nlceqx < -hugecon) = -hugecon;
            end
            if any(size([nlcineq; nlceq]) ~= size([nlcineqx; nlceqx])) || any(isnan([nlcineq; nlceq]) ~= isnan([nlcineqx; nlceqx])) || (~any(isnan([nlcineq; nlceq; nlcineqx; nlceqx])) && any(abs([0; nlcineq; nlceq] - [0; nlcineqx; nlceqx]) > cobyla_prec*max(1,abs([0; nlcineqx; nlceqx]))))
            % In the last few max of the above line, we put a 0 to avoid an empty result
                % Public/unexpected error
                error(sprintf('%s:InvalidConx', invoker), ...
                    '%s: UNEXPECTED ERROR: %s returns a con(x) that does not match x.', invoker, solver);
            end
        end

    end
end

% postpdfo ends
return
