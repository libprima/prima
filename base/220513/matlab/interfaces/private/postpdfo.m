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
% 1. All errors in this function are unexpected errors, which means they
% should not occur unless there is a bug in the code.
% 2. Some unexpected errors are public/external.
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
% N.B.: probinfo = [] if the invoker is a solver called by pdfo. In that case, postpdfo will return
% after verifying `output`, and then pdfo will call postpdfo again (after assigning `probinfo`) to
% do the postprocessing.
obligatory_probinfo_fields = {'raw_data', 'refined_data', 'fixedx', 'fixedx_value', ...
    'nofreex', 'infeasible_bound', 'infeasible_lineq', 'infeasible_leq', ...
    'trivial_lineq', 'trivial_leq', 'infeasible', 'scaled', 'scaling_factor', ...
    'shift', 'reduced', 'raw_type', 'raw_dim', 'refined_type', 'refined_dim', ...
    'feasibility_problem', 'user_options_fields', 'options', 'warnings', ...
    'hugenum', 'hugefun', 'hugecon'};
obligatory_options_fields = {'classical', 'debug', 'chkfunval', 'precision'};

% Who is calling this function? Is it a correct invoker?
invoker_list = ['pdfon', all_solvers()];
callstack = dbstack;
funname = callstack(1).name; % Name of the current function
if (length(callstack) == 1) || ~ismember(callstack(2).name, invoker_list)
    % Private/unexpected error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s should only be called by %s.', funname, funname, mystrjoin(invoker_list, ', '));
else
    invoker = callstack(2).name; % Name of the function who calls this function
end

% Verify the input before starting the real business
% Verify probinfo
if (length(callstack) >= 3) && strcmp(callstack(3).name, 'pdfon')
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

% Decide the solver that did the computation (needed for verifying output below).
if strcmp(invoker, 'pdfon')
    % In this case, the invoker is pdfo rather than a solver called by pdfo.
    % Thus probinfo is nonempty, and options has been read and verified as above.
    solver = options.solver;
else
    solver = invoker;
end
if isempty(solver) || ~ischarstr(solver) || ~ismember(solver, all_solvers())
    % Public/unexpected error
    error(sprintf('%s:InvalidSolver', invoker), '%s: UNEXPECTED ERROR: invalid solver passed to %s.', invoker, funname);
end

% Verify output
if ~isa(output, 'struct')
    % Public/unexpected error
    error(sprintf('%s:InvalidOutput', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an output that is not a structure', invoker, solver);
end
if ismember(solver, all_solvers('internal'))
    % For internal solvers, output should contain fhist, chist, and warnings
    obligatory_output_fields = [obligatory_output_fields, 'fhist', 'chist', 'warnings'];
end
if strcmp(solver, 'lincoan')
    % For lincoan, output should contain constr_modified
    obligatory_output_fields = [obligatory_output_fields, 'constr_modified'];
end
if ismember(solver, all_solvers('nonlinearly_constrained_solvers')) && ismember(solver, all_solvers('internal'))
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
% If the invoker is a solver called by pdfo, then let pdfo do the postprocessing.
% Put this after verifying output, because we will use the information in it.
if (length(callstack) >= 3) && strcmp(callstack(3).name, 'pdfon')
    x = output.x;
    fx = output.fx;
    exitflag = output.exitflag;
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% !!!!!! No reading of `probinfo` should happen before this line !!!!! %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% With the moderated extreme barrier (implemented when options.classical is false), all
% the function values that are NaN or larger than hugefun are replaced by hugefun;
% all the constraint values that are NaN or larger than hugecon are replaced by hugecon.
hugefun = probinfo.hugefun;
hugecon = probinfo.hugecon;

% Record solver name in output (do it after verifying that output is a structure).
output.algorithm = solver;

% Read information in output
x = output.x;  % x will be used later, for example, in length(x)
output = rmfield(output, 'x'); % output does not include x at return
fx = output.fx;
output = rmfield(output, 'fx'); % output does not include fx at return
exitflag = output.exitflag;
output = rmfield(output, 'exitflag'); % output does not include exitflag at return
nf = output.funcCount;
constrviolation = output.constrviolation;
if strcmp(solver, 'lincoan')
    constr_modified = output.constr_modified;
    output = rmfield(output, 'constr_modified');
end
if ~isfield(output, 'warnings') || isempty(output.warnings)
    output.warnings = {};
end

% Verify x
if ~isnumeric(x) || ~isreal(x) || ~isvector(x) || size(x,2) ~= 1
    % Public/unexpected error
    error(sprintf('%s:InvalidX', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an x that is not a real column or scalar.', invoker, solver);
end

% Verify fx
if ~isrealscalar(fx)
    % Public/unexpected error
    error(sprintf('%s:InvalidFx', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an fx that is not a real number.', invoker, solver);
end

% Verify exitflag
if ~isintegerscalar(exitflag)
    % Public/unexpected error
    error(sprintf('%s:InvalidExitFlag', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an exitflag that is not an integer', invoker, solver);
end

% Verify nf
if ~isintegerscalar(nf)
    % Public/unexpected error
    error(sprintf('%s:InvalidNF', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an nf that is not an integer.', invoker, solver);
end
if nf <= 0
    % If prepdfo works properly, then nf <= 0 should never happen.
    % Public/unexpected error
    error(sprintf('%s:InvalidNF', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns nf = 0 unexpectedly with exitflag %d.', invoker, solver, exitflag);
end

% For internal solvers:
% xhist is either empty or containing the last nhist iterates of the solver;
% nlcihist is either empty or containing the nonlinear inequality constraint values of the
% last nhist iterates of the solver;
% nlcehist is either empty or containing the nonlinear equality constraint values of the
% last nhist iterates of the solver;
% fhist contains the function values of the last nhist iterates of the solver.
if isfield(output, 'fhist')
    nhist = length(output.fhist);
else
    nhist = 0;
end

% Read and verify xhist
output_has_xhist = isfield(output, 'xhist');
if isfield(output, 'xhist')
    xhist = output.xhist;
else
    xhist = [];
end
if ~isempty(xhist) && (~isrealmatrix(xhist) || any(size(xhist) ~= [length(x), nhist]))  % x is set before
    % Public/unexpected error
    error(sprintf('%s:InvalidXhist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an xhist that is not a real matrix of size (n, min(nf, maxhist)).', invoker, solver);
end

% Read and verify fhist
if isfield(output, 'fhist')
    fhist = output.fhist;
else % External solvers may not return fhist
    fhist = [];
end
if ~isempty(fhist) && ~(isrealvector(fhist) && length(fhist) == nhist)
    % Public/unexpected error
    error(sprintf('%s:InvalidFhist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an fhist that is not a real vector of length min(nf, maxhist).', invoker, solver);
end
if ~options.classical && ~probinfo.infeasible && ~probinfo.nofreex
    if any(fhist > hugefun) || any(isnan(fhist))
        % Public/unexpected error
        error(sprintf('%s:InvalidFhist', invoker), ...
             '%s: UNEXPECTED ERROR: %s returns an fhist with NaN or values larger than hugefun = %g; this is impossible except in the classical mode.', invoker, solver, hugefun);
    elseif ~isempty(fhist) && max(fhist) == hugefun
        wid = sprintf('%s:ExtremeBarrier', invoker);
        wmsg = sprintf('%s: the moderated extreme barrier is invoked; function values that are NaN or larger than hugefun = %g are replaced by hugefun.', invoker, hugefun);
        warning(wid, '%s', wmsg);
        output.warnings = [output.warnings, wmsg];
    end
end

% If the problem is a feasibility problem, set fx to [], and remove fhist from output.
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
if ~isrealscalar(constrviolation)
    % Public/unexpected error
    error(sprintf('%s:InvalidConstrViolation', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns a constrviolation that is not a real number.', invoker, solver)
end

% Read and verify chist
output_has_chist = isfield(output, 'chist');
if isfield(output, 'chist')
    chist = output.chist;
else % External solvers may not return chist
    chist = constrviolation + zeros(1, nhist);
end
if ~(isempty(chist) && ismember(solver, all_solvers('without_constraints'))) && ~(isrealvector(chist) && length(chist) == nhist)
    % Public/unexpected error
    error(sprintf('%s:InvalidChist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns a chist that is not a real vector of length min(nf, maxfhist).', invoker, solver);
end
if ~options.classical && ~probinfo.infeasible && ~probinfo.nofreex
    if strcmp(solver, 'cobylan') && (any(chist > hugecon) || any(isnan(chist)))
        % Public/unexpected error
        error(sprintf('%s:InvalidChist', invoker), ...
             '%s: UNEXPECTED ERROR: %s returns a chist with NaN or values larger than hugecon = %g; this is impossible except in the classical mode.', invoker, solver, hugecon);
    elseif ~isempty(chist) && max(chist) == hugecon
        wid = sprintf('%s:ExtremeBarrier', invoker);
        wmsg = sprintf('%s: the moderated extreme barrier is invoked; constraint values that are NaN or larger than hugecon = %g are replaced by hugecon.', invoker, hugecon);
        warning(wid, '%s', wmsg);
        output.warnings = [output.warnings, wmsg];
    end
end

% Read and verify nlcineq and nlceq
if isfield(output, 'nlcineq')
    nlcineq = output.nlcineq;
else
    nlcineq = [];  % Needed later, e.g., to decide cstrv.
end
if isfield(output, 'nlceq')
    nlceq = output.nlceq;
else
    nlceq = [];  % Needed later, e.g., to decide cstrv.
end
if ~strcmp(probinfo.refined_type, 'nonlinearly-constrained') && ~(isempty(nlcineq) && isempty(nlceq))
    % Public/unexpected error
    error(sprintf('%s:InvalidNonlinearConstraint', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns values of nonlinear constraints for a problem without such constraints.', invoker, solver);
end
if (isfield(output, 'nlcineq') && ~isfield(output, 'nlceq')) || (~isfield(output, 'nlcineq') && isfield(output, 'nlceq'))
    % Public/unexpected error
    error(sprintf('%s:InvalidNonlinearConstraint', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns only one of nlcineq and nlceq; it should return both of them or neither of them.', invoker, solver);
end

% Read and verify nlcihist and nlcehist
if isfield(output, 'nlcihist')
    nlcihist = output.nlcihist;
else
    nlcihist = [];  % Must be defined.
end
if isfield(output, 'nlcehist')
    nlcehist = output.nlcehist;
else
    nlcehist = [];  % Must be defined.
end
if ~strcmp(probinfo.refined_type, 'nonlinearly-constrained') && ~(isempty(nlcihist) && isempty(nlcehist))
    % Public/unexpected error
    error(sprintf('%s:InvalidNonlinearConstraint', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns history of nonlinear constraints for a problem without such constraints.', invoker, solver);
end
if (isfield(output, 'nlcihist') && ~isfield(output, 'nlcehist')) || (~isfield(output, 'nlcihist') && isfield(output, 'nlcehist'))
    % Public/unexpected error
    error(sprintf('%s:InvalidNonlinearConstraint', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns only one of nlcihist and nlcehist; it should return both of them or neither of them.', invoker, solver);
end
if ~isempty(nlcihist) && (~isrealmatrix(nlcihist) || any(size(nlcihist) ~= [length(nlcineq), nhist]))
    % Public/unexpected error
    error(sprintf('%s:InvalidNlcihist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns a nlcihist that is not a real matrix of min(nf, maxfhist) rows or empty.', invoker, solver);
end
if ~isempty(nlcehist) && (~isrealmatrix(nlcehist) || any(size(nlcehist) ~= [length(nlceq), nhist]))
    % Public/unexpected error
    error(sprintf('%s:InvalidNlcehist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns a nlcehist that is not a real matrix of min(nf, maxfhist) rows or empty.', invoker, solver);
end


% After verification, extract and process the data.

% The problem was (possibly) scaled. Scale it back.
% The scaling affects constrviolation when there are bound constraint.
% Hence constrviolation has to be recalculated so that it equals the
% constraint violation of the returned x with respect to the original problem.
% Ideally, chist should also be recalculated. However, it is impossible
% because we do not save the history of x. Therefore, when
% probinfo.scaled is true, chist is not the history of constraint violation
% of the original problem but the scaled one. It it not consistent with
% constrviolation. Without saving of history of x, we cannot do better.
%
% Before recalculating constrviolation, save the one returned by the
% solver, because it will be used in debug mode when checking whether fx
% is consistent with fhist and chist. See the definition of fhistf for
% details.
cstrv_returned = constrviolation;
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

    % Scale x back
    x = probinfo.scaling_factor.*x + probinfo.shift;
    % Scale xhist back
    xhist = probinfo.scaling_factor.*xhist + probinfo.shift;

    % Scale bounds back
    lb = probinfo.scaling_factor.*probinfo.refined_data.lb + probinfo.shift;
    ub = probinfo.scaling_factor.*probinfo.refined_data.ub + probinfo.shift;

    % Calculate the constrviolation for the unscaled problem.
    constrviolation = get_cstrv(x, Aineq, bineq, Aeq, beq, lb, ub, nlcineq, nlceq);
end

% The problem was (possibly) reduced. Get the full x and xhist.
if probinfo.reduced
    freex_value = x;
    x = NaN(length(x)+length(probinfo.fixedx_value), 1);
    x(probinfo.fixedx) = probinfo.fixedx_value;
    x(~probinfo.fixedx) = freex_value;

    freexhist = xhist;
    xhist = NaN(length(x), size(xhist, 2));  % x has been recovered; size(xhist,2) may not be nhist but 0.
    xhist(probinfo.fixedx, :) = probinfo.fixedx_value*ones(1, size(xhist, 2));
    xhist(~probinfo.fixedx, :) = freexhist;
end

% Include the possibly revised xhist to output if necessary.if
if output_has_xhist
    output.xhist = xhist;
end

% Set output.constrviolation to the revised constraint violation.
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

% Record the return message in output.message according to exitflag
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
case 20
    output.message = sprintf('Return from %s because the trust-region iteration has been performed maxtr (= 2*maxfun) times.', invoker);
case -1
    output.message = sprintf('Return from %s because NaN occurs in x.', solver);
case -2  % This cannot happen if the moderated extreme barrier is implemented, which is the case when options.classical is false.
    if strcmp(solver, 'cobylan')
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
    if ismember(solver, all_solvers('with_constraints'))
        fhistf = fhistf(chist <= max(cstrv_returned, 0));
    end
    minf = min([fhistf, fx]);
%%%% Zaikun 2021-05-26: The following test is disabled for lincoa for the moment. lincoa may not pass it.
%%%%    if (fx ~= minf) && ~(isnan(fx) && isnan(minf)) && ~(strcmp(solver, 'lincoan') && constr_modified)
    if (fx ~= minf) && ~(isnan(fx) && isnan(minf)) && ~strcmp(solver, 'lincoan')
        % Public/unexpected error
        error(sprintf('%s:InvalidFhist', invoker), ...
             '%s: UNEXPECTED ERROR: %s returns an fhist that does not match nf or fx.', invoker, solver);
    end

    % Check whether constrviolation is correct
    cobylan_prec = 1e-6;
    lincoan_prec = 1e-9;
    % COBYLA cannot ensure fx == fun(x) or constr == con(x) due to rounding
    % errors. Instead of checking the equality, we check whether the
    % relative error is within cobylan_prec.
    % There can also be a difference between constrviolation and cstrv due
    % to rounding errors, especially if the problem is scaled.
    constrviolation = 0;
    if isfield(output, 'constrviolation')
        constrviolation = output.constrviolation;
    end
    if strcmp(solver, 'bobyqan') && (max([chist, constrviolation]) > 0) && ~probinfo.infeasible
        % Public/unexpected error
        error(sprintf('%s:InvalidChist', invoker), ...
             '%s: UNEXPECTED ERROR: %s is a feasible solver yet it returns positive constrviolations.', invoker, solver);
    end
    if strcmp(options.precision, 'double') && ((strcmp(solver, 'lincoan') && ~constr_modified) || strcmp(solver, 'cobylan'))
        Aineq = probinfo.raw_data.Aineq;
        bineq = probinfo.raw_data.bineq;
        Aeq = probinfo.raw_data.Aeq;
        beq = probinfo.raw_data.beq;
        lb = probinfo.raw_data.lb(:);
        ub = probinfo.raw_data.ub(:);
        cstrv = get_cstrv(x, Aineq, bineq, Aeq, beq, lb, ub, nlcineq, nlceq);
        if ~(isnan(cstrv) && isnan(constrviolation)) && ~(cstrv == inf && constrviolation == inf) && ~(abs(constrviolation-cstrv) <= lincoan_prec*max(1,abs(cstrv)) && strcmp(solver, 'lincoan')) && ~(abs(constrviolation-cstrv) <= cobylan_prec*max(1,abs(cstrv)) && strcmp(solver, 'cobylan'))
            % Public/unexpected error
            error(sprintf('%s:InvalidChist', invoker), ...
              '%s: UNEXPECTED ERROR: %s returns a constrviolation that does not match x.', invoker, solver);
        end

        if isnan(fx)
            cf = chist(isnan(fhist));
        else
            cf = chist(fhist == fx);
        end
        if (nhist >= nf) && ~any(cf == cstrv_returned) && ~(isnan(cstrv_returned) && ~any(~isnan(cf)))
            % Public/unexpected error
            % Note: When nhist < nf, FHIST and CHIST do not contain the whole history.
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
        % Due to the moderated extreme barrier (implemented when options.classical is false),
        % all function values that are NaN or larger than hugefun are replaced by hugefun.
        if (funx ~= funx) || (funx > hugefun)
            funx = hugefun;
        end
        %if (funx ~= fx) && ~(isnan(fx) && isnan(funx))
        % it seems that COBYLA can return fx ~= fun(x) due to rounding
        % errors. Therefore, we cannot use "fx ~= funx" to check COBYLA
        if ~(isnan(fx) && isnan(funx)) && ~((fx == funx) || (abs(funx-fx) <= cobylan_prec*max(1, abs(fx)) && strcmp(solver, 'cobylan')))
            % Public/unexpected error
            error(sprintf('%s:InvalidFx', invoker), ...
                '%s: UNEXPECTED ERROR: %s returns an fx that does not match x.', invoker, solver);
        end

        % Check whether fhist = fun(xhist)
        if ~isempty(fhist) && ~isempty(xhist)
            fhistx = zeros(1, nhist);  % When the objective is empty, the objective function used in computation was 0.
            if ~isempty(objective)
                for k = 1 : nhist
                    fhistx(k) = objective(xhist(:, k));
                end
            end
            % Due to the moderated extreme barrier (implemented when options.classical is false),
            % all function values that are NaN or above hugefun are replaced by hugefun.
            fhistx(fhistx ~= fhistx | fhistx > hugefun) = hugefun;
            if any(~(isnan(fhist) & isnan(fhistx)) & ~((fhist == fhistx) | (abs(fhistx-fhist) <= cobylan_prec*max(1, abs(fhist)) & strcmp(solver, 'cobylan'))))
                % Public/unexpected error
                error(sprintf('%s:InvalidFx', invoker), ...
                    '%s: UNEXPECTED ERROR: %s returns an fhist that does not match xhist.', invoker, solver);
            end
        end

        % Check whether [output.nlcineq,  output.nlceq] = nonlcon(x)
        nlcineqx = [];
        nlceqx = [];
        nonlcon = probinfo.raw_data.nonlcon;
        if ~isempty(nonlcon)
            [nlcineqx, nlceqx] = feval(nonlcon, x);
            % Due to the moderated extreme barrier (implemented when options.classical is false),
            % all constraint values that are NaN or above hugecon are replaced by hugecon.
            nlcineqx(nlcineqx ~= nlcineqx | nlcineqx > hugecon) = hugecon;
            % All constraint values below -hugecon are replaced by -hugecon to avoid numerical difficulties.
            nlcineqx(nlcineqx < -hugecon) = -hugecon;
            nlceqx(nlceqx ~= nlceqx | nlceqx > hugecon) = hugecon;
            nlceqx(nlceqx < -hugecon) = -hugecon;
        end
        if any(size([nlcineq; nlceq]) ~= size([nlcineqx; nlceqx])) || any(isnan([nlcineq; nlceq]) ~= isnan([nlcineqx; nlceqx])) || (~any(isnan([nlcineq; nlceq; nlcineqx; nlceqx])) && any(abs([0; nlcineq; nlceq] - [0; nlcineqx; nlceqx]) > cobylan_prec*max(1,abs([0; nlcineqx; nlceqx]))))
        % In the last few max of the above line, we put a 0 to avoid an empty result
            % Public/unexpected error
            error(sprintf('%s:InvalidConx', invoker), ...
                '%s: UNEXPECTED ERROR: %s returns a con(x) that does not match x.', invoker, solver);
        end

        if ~isempty(xhist) && (isfield(output, 'nlcihist') || isfield(output, 'nlcehist') || isfield(output, 'chist'))
            % Calculate nlcihistx and nlcehistx according to xhist.
            nlcihistx = [];
            nlcehistx = [];
            nonlcon = probinfo.raw_data.nonlcon;
            if ~isempty(nonlcon)
                nlcihistx = NaN(length(nlcineq), nhist);
                nlcehistx = NaN(length(nlceq), nhist);
                for k = 1 : nhist
                    [nlcihistx(:, k), nlcehistx(:, k)] = feval(nonlcon, xhist(:, k));
                end
                % Due to the moderated extreme barrier (implemented when options.classical is false),
                % all constraint values that are NaN or above hugecon are replaced by hugecon.
                nlcihistx(nlcihistx ~= nlcihistx | nlcihistx > hugecon) = hugecon;
                % All constraint values below -hugecon are replaced by -hugecon to avoid numerical difficulties.
                nlcihistx(nlcihistx < -hugecon) = -hugecon;
                nlcehistx(nlcehistx ~= nlcehistx | nlcehistx > hugecon) = hugecon;
                nlcehistx(nlcehistx < -hugecon) = -hugecon;
            end

            % Check whether [output.nlcihist, output.nlcehist] = nonlcon(xhist).
            % Do not replace the condition by isempty([nlcihist; nlcehist]), as it may hold when it should not!
            if isfield(output, 'nlcihist') || isfield(output, 'nlcehist')
                if any(size([nlcihist; nlcehist]) ~= size([nlcihistx; nlcehistx])) || ...
                        any(isnan([nlcihist; nlcehist]) ~= isnan([nlcihistx; nlcehistx]), 'all') || ...
                        (~any(isnan([nlcihist; nlcehist; nlcihistx; nlcehistx]), 'all') && ...
                        any(abs([zeros(1, nhist); nlcihist; nlcehist] - [zeros(1, nhist); nlcihistx; nlcehistx]) ...
                        > cobylan_prec*max(1,abs([zeros(1, nhist); nlcihistx; nlcehistx])), 'all'))
                    % In the last few max of the above line, we put a 0 to avoid an empty result
                    % Public/unexpected error
                    error(sprintf('%s:InvalidConx', invoker), ...
                        '%s: UNEXPECTED ERROR: %s returns an nlcihist or nlcehist that does not match xhist.', invoker, solver);
                end
            end

            % Check whether chist = constrviolation(xhist).
            if (strcmp(solver, 'lincoan') && ~constr_modified) || strcmp(solver, 'cobylan')
                % In this case, chist has been set to output.chist, but output.chist has been
                % removed if the problem is unconstrained.
                chistx = NaN(1, nhist);
                for k = 1 : nhist
                    if isempty(nonlcon)
                        chistx(k) = get_cstrv(xhist(:, k), Aineq, bineq, Aeq, beq, lb, ub, [], []);
                    else
                        chistx(k) = get_cstrv(xhist(:, k), Aineq, bineq, Aeq, beq, lb, ub, nlcihistx(:, k), nlcehistx(:, k));
                    end
                end
                if any(~(isnan(chist) & isnan(chistx)) & ~((chist == chistx) | abs(chistx-chist) <= lincoan_prec*max(1, abs(chist)) & strcmp(solver, 'lincoan') | (abs(chistx-chist) <= cobylan_prec*max(1, abs(chist)) & strcmp(solver, 'cobylan'))))
                    % Public/unexpected error
                    error(sprintf('%s:InvalidFx', invoker), ...
                        '%s: UNEXPECTED ERROR: %s returns an chist that does not match xhist.', invoker, solver);
                end
            end
        end
    end % chkfunval ends
end

% postpdfo ends
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cstrv = get_cstrv(x, Aineq, bineq, Aeq, beq, lb, ub, nlcineq, nlceq)
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
cstrv = max([0; rineq; abs(req); lb-x; x-ub; nlcineq; abs(nlceq)], [], 'includenan');
% max(X, [], 'includenan') returns NaN if X contains NaN, and maximum of X otherwise
return
