function [x, fx, exitflag, output] = postpdfo(x, fx, exitflag, output, nf, fhist, constrviolation, chist, options, probinfo)
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

% With extreme barrier (implemented when options.classical=false), all
% the function values that are NaN or larger than hugefun are replaced
% by hugefun; all the constraint values that are NaN or larger than
% hugecon are replaced by hugecon. hugefun and hugecon are defined in
% pdfoconst.F, and can be obtained by gethuge. 
hugefun = gethuge('fun'); 
hugecon = gethuge('con');

% Who is calling this function? Is it a correct invoker?
invoker_list = {'pdfo', 'uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla'};
unconstrained_solver_list = {'uobyqa', 'newuoa'};
constrained_solver_list = {'bobyqa', 'lincoa', 'cobyla'};
solver = options.solver; % The solver that did the computation
callstack = dbstack;
funname = callstack(1).name; % Name of the current function 
if (length(callstack) == 1) || ~ismember(callstack(2).name, invoker_list) 
    % Private/unexpcted error
    error(sprintf('%s:InvalidInvoker', funname), ...
    '%s: UNEXPECTED ERROR: %s should only be called by %s.', funname, funname, mystrjoin(invoker_list, ', '));
else
    invoker = callstack(2).name; % Name of the function who calls this function
end

% If the invoker is a solver called by pdfo, then let pdfo do the postprecessing
if (length(callstack) >= 3) && strcmp(callstack(3).name, 'pdfo')
    output.funcCount = nf;
    output.constrviolation = constrviolation;
    output.fhist = fhist;
    output.chist = chist;
    return
end

% Varify the input before starting the real business
if ~isnumeric(x) || ~isreal(x) || ~isvector(x) || size(x,2)~=1
    % Public/unexpcted error
    error(sprintf('%s:InvalidX', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an x that is not a real column or scalar.', invoker, solver);
end

if ~isnumeric(fx) || ~isreal(fx) || ~isscalar(fx)
    % Public/unexpcted error
    error(sprintf('%s:InvalidFx', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an fx that is not a real number.', invoker, solver);
end

if ~isnumeric(exitflag) || ~isscalar(exitflag) || ~isreal(exitflag) || rem(exitflag, 1)~=0
    % Public/unexpcted error
    error(sprintf('%s:InvalidExitFlag', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an exitflag that is not an integer', invoker, solver);
end

if ~isa(output, 'struct')
    % Public/unexpcted error
    error(sprintf('%s:InvalidOutput', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an output that is not a structure', invoker, solver);
end

if ~isnumeric(nf) || ~isscalar(nf) || ~isreal(nf) || rem(nf, 1)~=0 || nf < 0 
    % Public/unexpcted error
    error(sprintf('%s:InvalidNF', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an nf that is not a nonnegative integer.', invoker, solver);
end

if (nf == 0) && ~probinfo.infeasible 
    % If prepdfo works properly, then nf=0 should not happen unless infeasibility is detecte 
    % Public/unexpcted error
    error(sprintf('%s:InvalidNF', invoker), ...
    '%s: UNEXPECTED ERROR: %s returns nf=0 unexpectedly.', invoker, solver);
end

if ~isempty(fhist) && (~isnumeric(fhist) || ~isreal(fhist) || ~isvector(fhist) || (nf ~= length(fhist)))
    % Public/unexpected error
    error(sprintf('%s:InvalidFhist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an fhist that is not a real vector of length nf.', invoker, solver);
end
if ~options.classical
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

if ~isnumeric(constrviolation) || ~isscalar(constrviolation) || ~isreal(constrviolation)
    % Public/unexpected error
    error(sprintf('%s:InvalidConstrViolation', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns a constrviolation that is not a real number.', invoker, solver)
end

if ~(isempty(chist) && (ismember(solver, unconstrained_solver_list) || nf == 0)) && (~isnumeric(chist) || ~isreal(chist) || ~isvector(chist) || (nf ~= length(chist)))
    % Public/unexpected error
    error(sprintf('%s:InvalidChist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns a chist that is not a real vector of length nf.', invoker, solver);
end
if ~options.classical
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

if ~isa(options, 'struct')
    % Public/unexpected error
    error(sprintf('%s:InvalidChist', invoker), ...
        '%s: UNEXPECTED ERROR: %s returns an options that is not a structure.', invoker, solver);
end

if isfield(probinfo, 'warnings') && ~isempty(probinfo.warnings) && ~(isa(probinfo.warnings, 'cell') && isvector(probinfo.warnings))
    % Public/unexpected error
    error(sprintf('%s:InvalidWarnings', invoker),...
    '%s: UNEXPECTED ERROR: %s returns warnings that is not a cell of strings.', invoker, solver);
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
        constrviolation = max(max(0, [rineq; req; -req; lb-x; x-ub])./max(1, abs([bineq; beq; -beq; lb; ub]))); 
        % LINCOA returns a relative constraint violation
    else
        nlcineq = [];
        nlceq = [];
        if isfield(output, 'nlcineq') || isfield(output, 'nlceq')
            nlcineq = output.nlcineq;
            nlceq = output.nlceq;
        end
        constrviolation = max([0; rineq; req; -req; lb-x; x-ub; nlcineq; nlceq; -nlceq]);
        if any(isnan([rineq; req; -req; lb-x; x-ub; nlcineq; nlceq; -nlceq]))
            constrviolation = NaN;
        end
        % COBYLA returns an absolute constraint violation (there is
        % nothing to compare with, because con is a blackbox)
        % In the max, we put a 0 to avoid an empty result
    end
end

% The problem was (possibly) reduced. Get the full x.
if probinfo.reduced 
    freex_value = x;
    x = NaN(length(x)+length(probinfo.fixedx_value), 1);
    x(probinfo.fixedx) = probinfo.fixedx_value;
    x(~probinfo.fixedx) = freex_value;
end

% Set output.{nf, constrviolation, fhist, chist, solver}
output.constrviolation = constrviolation;
output.funcCount = nf;
output.fhist = fhist;
output.chist = chist;
output.algorithm = solver;

% Revise constrviolation and chist according to problem type
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
elseif nf > 0 && strcmp(probinfo.refined_type, 'unconstrained') && ~strcmp(probinfo.raw_type, 'unconstrained') 
    output.constrviolation = 0.0;
    output.chist = zeros(1, nf);
end

% Record the retrun message in output.message according to exitflag 
switch exitflag % If prepdfo works properly, then 5, 6, 10, 11, 12 should never happen
case 0
    output.message = sprintf('Return from %s because the lower bound for the trust region radius is reached.', solver);
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

% Record the warnings in output.warnings
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

% More careful checks about fx, constrviolation, fhist, and chist.
% We do this only if the code is in debug mode but not in classical
% mode. The classical mode cannot pass these checks.  
if options.debug && ~options.classical && (nf > 0) 
    % Check whether fx is 'optimal' 
    fhistf = fhist;
    if ismember(solver, constrained_solver_list) 
        fhistf = fhistf(chist <= max(constrv_returned, 0));
    end
    minf = min([fhistf,fx]);
    if (fx ~= minf) && ~(isnan(fx) && isnan(minf)) && ~(strcmp(solver, 'lincoa') && output.constr_modified)
        % Public/unexpected error
        error(sprintf('%s:InvalidFhist', invoker), ...
             '%s: UNEXPECTED ERROR: %s returns an fhist that does not match nf or fx.', invoker, solver);
    end

    % Check whether constrviolation is correct
    cobyla_prec = 1e-12; 
    lincoa_prec = 1e-14; 
    % COBYLA cannot ensure fx=fun(x) or conval=con(x) due to rounding
    % errors. Instead of checking the equality, we check whether the
    % relative error is within cobyla_prec. 
    constrviolation = 0;
    if isfield(output, 'constrviolation')
        constrviolation = output.constrviolation;
    end
    if strcmp(solver, 'bobyqa') && (max([chist, constrviolation]) > 0) 
        % Public/unexpected error
        error(sprintf('%s:InvalidChist', invoker), ...
             '%s: UNEXPECTED ERROR: %s is a feasible solver yet it returns positive constrviolations.', invoker, solver);
    end
    if (strcmp(solver, 'lincoa') && ~output.constr_modified) || strcmp(solver, 'cobyla')
        Aineq = probinfo.raw_data.Aineq;
        bineq = probinfo.raw_data.bineq;
        Aeq = probinfo.raw_data.Aeq;
        beq = probinfo.raw_data.beq;
        lb = probinfo.raw_data.lb;
        ub = probinfo.raw_data.ub;
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
            conv = max(max(0, [rineq; req; -req; lb-x; x-ub])./max(1, abs([bineq; beq; -beq; lb; ub]))); 
            % LINCOA returns a relative constraint violation
        else
            nlcineq = [];
            nlceq = [];
            if isfield(output, 'nlcineq') || isfield(output, 'nlceq')
                nlcineq = output.nlcineq;
                nlceq = output.nlceq;
            end
            conv = max([0; rineq; req; -req; lb-x; x-ub; nlcineq; nlceq; -nlceq]);
            if any(isnan([rineq; req; -req; lb-x; x-ub; nlcineq; nlceq; -nlceq]))
                conv = NaN;
            end
            % COBYLA returns an absolute constraint violation (there is
            % nothing to compare with, because con is a blackbox)
            % In the max, we put a 0 to avoid an empty result
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
        nlcineq = [];
        nlceq = [];
        if (isfield(output, 'nlcineq'))
            nlcineq = output.nlcineq;
        end
        if (isfield(output, 'nlceq'))
            nlceq = output.nlceq;
        end
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

if strcmp(solver, 'lincoa')
    output = rmfield(output, 'constr_modified');
end

% Print the return message if options.quiet=false
if ~options.quiet
    fprintf('\n%s\n\n', output.message)
end

warning('on', 'backtrace'); % Recover the default warning behavior of displaying stack trace, which was disabled by pdfo or its solvers

% postpdfo ends
return
