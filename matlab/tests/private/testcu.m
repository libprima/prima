function [mrec, mmin, output] = testcu(solvers, options)

% Default options
rhobeg = 1;
rhoend = 1e-6;
maxfun_dim = 500;
%maxfun = 50000;
maxfun = 20000;
maxit = 1000;
ftarget = -inf;
perm = false;
randomizex0 = 0;
eval_options = struct();
nr = 5;
ctol = 1e-10;
cpenalty = 1e10;
%ctol = inf
%cpenalty  = 0
type = 'ubln'; % The types of problems to test
mindim = 1; % The minimal dimension of problems to test
solvers = lower(solvers);
if startsWith(solvers{1}, 'cobyla') || startsWith(solvers{2}, 'cobyla')
    maxdim = 20; % The maximal dimension of problems to test
else
    maxdim = 50; % The maximal dimension of problems to test
end
mincon = 0; % The minimal number of constraints of problems to test
maxcon = min(5000, 100*maxdim); % The maximal number of constraints of problems to test
sequential = false;
%debug = false;
debug = true;
chkfunval = true;
output_xhist = true;
output_nlchist = true;
thorough_test = 0;
minip = 1;

% Set options
options = setopt(options, rhobeg, rhoend, maxfun_dim, maxfun, maxit, ftarget, perm, randomizex0, ...
    eval_options, nr, ctol, cpenalty, type, mindim, maxdim, mincon, maxcon, ...
    sequential, debug, chkfunval, output_xhist, output_nlchist, thorough_test, minip);

% Select the problems to test.
if isfield(options, 'list')
    plist = options.list;  % Only test problems in this list
    if (ischstr(plist))  % In case plist is indeed the name of a problem
        plist = {plist};
    end
else
    requirements = struct();
    requirements.mindim = options.mindim;
    requirements.maxdim = options.maxdim;
    requirements.mincon = options.mincon;
    requirements.maxcon = options.maxcon;
    requirements.type = options.type;
    requirements.blacklist = {};
    if startsWith(solvers{1}, 'cobyla') || startsWith(solvers{2}, 'cobyla')
        requirements.blacklist = [requirements.blacklist, {'CHEBYQADNE','HAIFAM','HIMMELBI','HYDCAR20','LUKSAN12','LUKSAN13','MSS1','SPANHYD','VANDERM1','VANDERM2','VANDERM3', 'TAX13322', 'TAXR13322'}]; % Takes more than 2 min to solve
        requirements.blacklist = [requirements.blacklist, {'DMN15102', 'DMN15103', 'DMN15332', 'DMN15333', 'DMN37142', 'DMN37143'}]; % Time-consuming
        requirements.blacklist = [requirements.blacklist, {'GMNCASE2'}];
        requirements.blacklist = [requirements.blacklist, {'VANDERM4', 'LAKES'}]; % The classical COBYLA encounters SIGFAULT
        requirements.blacklist = [requirements.blacklist, {'DALLASS', 'QPCBLEND'}]; % The profiling script on GitHub Actions seems to be blocked by them
    end
    if startsWith(solvers{1}, 'lincoa') || startsWith(solvers{2}, 'lincoa')
        requirements.blacklist = [requirements.blacklist, {'3PK', 'LSNNODOC', 'SIPOW3', 'SIPOW4', 'OET1', 'MAKELA4','TFI2', 'QPCBOEI2', 'QPNBOEI2'}]; % The classical LINCOA encounters SIGFAULT
        requirements.blacklist = [requirements.blacklist, {'TARGUS', 'ARGTRIGLS', 'VARDIM'}]; % Takes too long time
    end
    if startsWith(solvers{1}, 'uobyqa') || startsWith(solvers{2}, 'uobyqa')
        requirements.blacklist = [requirements.blacklist, {'BA-L1LS', 'BA-L1SPLS', 'CHNROSNB', 'CHNRSNBM', 'ERRINROS', 'ERRINRSM', 'TOINTGOR', 'TOINTPSP', 'VAREIGVL'}]; % Takes too long time
    end
    plist = secup(requirements);
end

%plist={'ALLINITU'};
np = length(plist);
ns = length(solvers);
nr = options.nr;
maxfun = options.maxfun;
sequential = options.sequential;
minip = options.minip;

% These arrays will record the function values and constraint values during the tests.
pdim = NaN(np, 1);  % Data profile needs the dimension of the problem.
frec = NaN(np, ns, nr, maxfun);
crec = NaN(np, ns, nr, maxfun);

% These arrays will record the reference function values and constraint values when there is an
% eval_options or `randomizex0` is positive.
fref = NaN(np, ns, maxfun);
cref = NaN(np, ns, maxfun);

permuted = options.perm;
has_eval_options = ~isempty(fieldnames(options.eval_options));
eval_options = options.eval_options;
randomizex0 = abs(options.randomizex0);
ref_options = rmfield(options, {'perm', 'randomizex0', 'eval_options'});

% `eval_options` and `randomizex0` can occur at the same time, but neither of them are compatible
% with `perm`.
assert(~permuted || ~(has_eval_options || randomizex0));


fprintf('\n\nThe testing options:\n')
display(options);
fprintf('\n\nThe evaluation options:\n')
display(eval_options);
if isfield(eval_options, 'noise')
    display(eval_options.noise);
end
if isfield(eval_options, 'dnoise')
    display(eval_options.dnoise);
end


if sequential
    for ip = minip : np
        orig_warning_state = warnoff(solvers);

        pname = plist{ip};

        fprintf('\n%3d. \t%s:\n', ip, upper(pname));

        prob = macup(pname);
        orig_prob = prob;
        prob.orig_objective = prob.objective;
        prob.orig_nonlcon = prob.nonlcon;
        prob.orig_x0 = prob.x0;
        pdim(ip) = length(prob.x0);

        if has_eval_options || randomizex0 > 0
            %fprintf('\nCalculate fref and cref\n');
            for is = 1 : ns
                [fref(ip, is, :), cref(ip, is, :)] = testsolv(solvers{is}, prob, ref_options);
            end
        end

        rng(ip); permutations = get_perms(nr, length(prob.x0));

        for ir = 1 : nr
            if has_eval_options
                prob.objective = @(x) evalfun(prob.orig_objective, x, eval_options, ir);
                if ~isempty(prob.orig_nonlcon)
                    prob.nonlcon = @(x) evalcon(prob.orig_nonlcon, x, eval_options, ir);
                end
            end

            if randomizex0 > 0
                rng(ir); r = randn(length(prob.x0), 1);
                prob.x0 = prob.orig_x0 + randomizex0*norm(prob.orig_x0)*r/norm(r);
            end

            if permuted
                prob = permprob(orig_prob, permutations(ir, :));
            end

            for is = 1 : ns
                [frec(ip, is, ir, :), crec(ip, is, ir, :)] = testsolv(solvers{is}, prob, options);
            end
        end

        warning(orig_warning_state); % Restore the behavior of displaying warnings
    end
else
    parfor ip = minip : np
        orig_warning_state = warnoff(solvers);

        pname = plist{ip};

        fprintf('\n%3d. \t%s:\n', ip, upper(pname));

        prob = macup(pname);
        orig_prob = prob;
        prob.orig_objective = prob.objective;
        prob.orig_nonlcon = prob.nonlcon;
        prob.orig_x0 = prob.x0;
        pdim(ip) = length(prob.x0);

        if has_eval_options || randomizex0 > 0
            %fprintf('\nCalculate fref and cref\n');
            for is = 1 : ns
                [fref(ip, is, :), cref(ip, is, :)] = testsolv(solvers{is}, prob, ref_options);
            end
        end

        rng(ip); permutations = get_perms(nr, length(prob.x0));

        for ir = 1 : nr
            if has_eval_options
                prob.objective = @(x) evalfun(prob.orig_objective, x, eval_options, ir);
                if ~isempty(prob.orig_nonlcon)
                    prob.nonlcon = @(x) evalcon(prob.orig_nonlcon, x, eval_options, ir);
                end
            end

            if randomizex0 > 0
                rng(ir); r = randn(length(prob.x0), 1);
                prob.x0 = prob.orig_x0 + randomizex0*norm(prob.orig_x0)*r/norm(r);
            end

            if permuted
                prob = permprob(orig_prob, permutations(ir, :));
            end

            for is = 1 : ns
                [frec(ip, is, ir, :), crec(ip, is, ir, :)] = testsolv(solvers{is}, prob, options);
            end
        end

        warning(orig_warning_state); % Restore the behavior of displaying warnings
    end
end


mrec = frec + options.cpenalty*crec;
mrec(crec > options.ctol) = NaN;
mrec(:,:,:,1) = frec(:,:,:,1) + options.cpenalty*crec(:,:,:,1); % Prevent mrec(:,:,:,1) from being NaN
mrec_min = min(min(min(mrec, [], 4), [], 3), [], 2);

if has_eval_options || randomizex0
    mref = fref + options.cpenalty*cref;
    mref_min = min(min(mref, [], 3), [], 2);
    mmin = min(mrec_min, mref_min);
else
    mmin = mrec_min;
end

output = struct();
output.plist = plist;
output.pdim = pdim;

return


function [fval_history, cv_history, output] = testsolv(solver, prob, options)

prob.options = setsolvopt(solver, length(prob.x0), options); % Set the options for the solver

if ischstr(solver)
    prob.options.classical = endsWith(solver, '_classical');
    if endsWith(solver, '_single')
        prob.options.precision = 'single';
    end
    if endsWith(solver, '_quadruple')
        prob.options.precision = 'quadruple';
    end
    % `regexprep` removes '_classical' in case 'solver' ends with it. Similar for '_single', '_quadruple'.
    solver = regexprep(solver, '_classical$', '');
    solver = regexprep(solver, '_single$', '');
    solver = regexprep(solver, '_quadruple$', '');
    solver = str2func(solver);
end

maxfun = options.maxfun;
fval_history = NaN(1, maxfun);
cv_history = NaN(1, maxfun);

has_eval_options = isfield(options, 'eval_options') && isstruct(options.eval_options) && ~isempty(fieldnames(options.eval_options));
prob.options.output_xhist = has_eval_options;

[~, ~, ~, output] = solver(prob);
% Some solvers (e.g., fmincon) may not respect maxfun. Indeed, PDFO solvers may also increase maxfun
% if it is too small (e.g., <= npt for NEWUOA).
nf = min(maxfun, output.funcCount);

if (nf >= 1)
    if has_eval_options
        xhist_cell = num2cell(output.xhist(:, 1:nf), 1);
        fval_history(1:nf) = cellfun(prob.orig_objective, xhist_cell);
        orig_cstrv = @(x) get_cstrv(x, prob.Aineq, prob.bineq, prob.Aeq, prob.beq, prob.lb, prob.ub, prob.orig_nonlcon);
        cv_history(1:nf) = cellfun(orig_cstrv, xhist_cell);
    else
        fval_history(1:nf) = output.fhist(1:nf);
        if isfield(output, 'chist')
            cv_history(1:nf) = max(0, output.chist(1:nf));
        else
            cv_history(1:nf) = zeros(1, nf);
        end
    end

    fval_history(nf+1:maxfun) = fval_history(nf);
    cv_history(nf+1:maxfun) = cv_history(nf);
else
    % Sometimes pdfo may return nf = 0, e.g., when it detects infeasibility.
    fval_history = prob.f0;
    cv_history = prob.constrv0;
end

return


function options = setopt(options, rhobeg, rhoend, maxfun_dim, maxfun, maxit, ftarget, perm, ...
        randomizex0, eval_options, nr, ctol, cpenalty, type, mindim, maxdim, mincon, maxcon, ...
        sequential, debug, chkfunval, output_xhist, output_nlchist, thorough_test, minip) % Set options

if (~isfield(options, 'rhoend'))
    options.rhoend = rhoend;
end
if (~isfield(options, 'rhobeg'))
    options.rhobeg = rhobeg;
end
if (~isfield(options, 'maxit'))
    options.maxit = maxit;
end
if (~isfield(options, 'ftarget'))
    options.ftarget = ftarget;
end
if (~isfield(options, 'ctol'))
    options.ctol = ctol;
end
if (~isfield(options, 'cpenalty'))
    options.cpenalty = cpenalty;
end
if ~isfield(options, 'perm')
    options.perm = perm;
end
options.perm = logical(options.perm);
if (~isfield(options, 'randomizex0'))
    options.randomizex0 = randomizex0;
end
options.randomizex0 = abs(options.randomizex0);
if (~isfield(options, 'nr'))
    options.nr = nr;
end
if (~isfield(options, 'type'))
    options.type = type;
end
if (~isfield(options, 'mindim'))
    options.mindim = mindim;
end
if (~isfield(options, 'maxdim'))
    options.maxdim = maxdim;
end
if (~isfield(options, 'mincon'))
    options.mincon = mincon;
end
if (~isfield(options, 'maxcon'))
    options.maxcon = maxcon;
end
if (~isfield(options, 'maxfun_dim'))
    options.maxfun_dim = maxfun_dim;
end
if (~isfield(options, 'maxfun'))
    options.maxfun = maxfun;
end
options.maxfun = min(options.maxfun, options.maxfun_dim*options.maxdim);
if (~isfield(options, 'sequential'))
    options.sequential = sequential;
end
if (~isfield(options, 'debug'))
    options.debug = debug;
end
if (~isfield(options, 'chkfunval'))
    options.chkfunval = chkfunval;
end
if (~isfield(options, 'output_xhist'))
    options.output_xhist = output_xhist;
end
if (~isfield(options, 'output_nlchist'))
    options.output_nlchist = output_nlchist;
end
if (~isfield(options, 'thorough_test'))
    options.thorough_test = thorough_test;
end
if (~isfield(options, 'minip'))
    options.minip = minip;
end

% Set eval_options
has_eval_options = isfield(options, 'eval_options') && isstruct(options.eval_options) && ~isempty(fieldnames(options.eval_options));
if ~has_eval_options
    options.eval_options = eval_options;
end

if has_eval_options
    eval_options = options.eval_options;

    noise.type = 'relative';
    noise.nature = 'normal';
    noise.level = 0;
    if isfield(eval_options, 'noise') && isnumeric(eval_options.noise) && isscalar(eval_options.noise)
        noise.level = abs(eval_options.noise);
    elseif isfield(eval_options, 'noise') && isstruct(eval_options.noise)
        noise = eval_options.noise;
        if ~isfield(noise, 'type')
            noise.type = 'relative';
        end
        if ~isfield(noise, 'nature')
            noise.nature = 'normal';
        end
        if ~isfield(noise, 'level')
            noise.level = 1e-3;  % The default noise level if `noise` is present in `eval_options`
        end
        noise.level = abs(noise.level);
    end
    eval_options.noise = noise;
    if eval_options.noise.level == 0
        eval_options = rmfield(eval_options, 'noise');
    end

    dnoise.type = 'relative';
    dnoise.level = 0;
    if isfield(eval_options, 'dnoise') && isnumeric(eval_options.dnoise) && isscalar(eval_options.dnoise)
        dnoise.level = abs(eval_options.dnoise);
    elseif isfield(eval_options, 'dnoise') && isstruct(eval_options.dnoise)
        dnoise = eval_options.dnoise;
        if ~isfield(dnoise, 'type')
            dnoise.type = 'relative';
        end
        if ~isfield(dnoise, 'level')
            dnoise.level = 1e-3;  % The default dnoise level if `dnoise` is present in `eval_options`
        end
        dnoise.level = abs(dnoise.level);
    end
    eval_options.dnoise = dnoise;
    if eval_options.dnoise.level == 0
        eval_options = rmfield(eval_options, 'dnoise');
    end

    if isfield(eval_options, 'signif1')
        eval_options.signif = 1;
    elseif isfield(eval_options, 'signif2')
        eval_options.signif = 2;
    elseif isfield(eval_options, 'signif3')
        eval_options.signif = 3;
    elseif isfield(eval_options, 'signif4')
        eval_options.signif = 4;
    elseif isfield(eval_options, 'signif5')
        eval_options.signif = 5;
    elseif isfield(eval_options, 'signif6')
        eval_options.signif = 6;
    end

    if isfield(eval_options, 'single')
        eval_options.single = true;
    end

    options.eval_options = eval_options;
end

eval_options = options.eval_options;

% Revise options.nr
noisy_eval = (isfield(eval_options, 'noise') && eval_options.noise.level > 0);
if ~(options.perm || options.randomizex0 > 0 || noisy_eval)
    options.nr = 1;
end

% Revise options.ctol and options.cpenalty
if isfield(eval_options, 'dnoise')
    options.ctol = max(options.ctol, eval_options.dnoise.level);
    options.cpenalty = min(options.cpenalty, 100/options.ctol);
end
if isfield(eval_options, 'noise')
    options.ctol = max(options.ctol, eval_options.noise.level);
    options.cpenalty = min(options.cpenalty, 100/options.ctol);
end
if isfield(eval_options, 'signif')
    options.ctol = max(options.ctol, 10^(-eval_options.signif));
    options.cpenalty = min(options.cpenalty, 100/options.ctol);
end
if isfield(eval_options, 'single') && eval_options.single
    options.ctol = max(options.ctol, 1e-5);
    options.cpenalty = min(options.cpenalty, 100/options.ctol);
end
if options.randomizex0 > 0
    options.ctol = max(options.ctol, 1e-5);
    options.cpenalty = min(options.cpenalty, 100/options.ctol);
end

return



function solv_options = setsolvopt(solv, n, options)

solv_options = struct();
solv_options.rhobeg = options.rhobeg;
solv_options.rhoend = options.rhoend;
solv_options.maxfun = min(options.maxfun_dim*n, options.maxfun);
solv_options.ftarget = options.ftarget;
solv_options.output_xhist = options.output_xhist;
solv_options.output_nlchist = options.output_nlchist;
solv_options.iprint = 0;
solv_options.quiet = true;
solv_options.debug = options.debug;
solv_options.chkfunval = options.chkfunval;
%solv_options.scale = true;

if (strcmpi(solv, 'fmincon'))
    solv_options = optimoptions('fmincon');
    solv_options.MaxFunctionEvaluations = min(options.maxfun_dim*n, options.maxfun);
    solv_options.MaxIterations = options.maxit;
    solv_options.ObjectiveLimit = options.ftarget;
    solv_options.OptimalityTolerance = options.rhoend;
    solv_options.StepTolerance = options.rhoend;
    solv_options.ConstraintTolerance = min(1e-6, options.rhoend);
end
return


function f = evalf(f, x, options)

if isfield(options, 'noise')
    noise = options.noise;
    if isstruct(noise) && isfield(noise, 'level') && noise.level > 0
        seed = 0.3*sin(1e8*abs(f))+0.3*cos(1e8*norm(x,9)) + 0.3*sin(100*norm(x,1))*cos(100*norm(x,Inf)) + 0.1*cos(norm(x));
        rng(min(options.ir*ceil(abs(10e6*seed)), 2^31));  % rng accepts integers between 0 and 2^32 - 1.

        switch lower(noise.nature)
        case {'uniform', 'u'}
            r = 2*rand-1;
        otherwise
            r = randn;
        end

        switch lower(noise.type)
        case {'absolute', 'additive', 'add', 'a', '+'}
            f = f + noise.level*r;
        otherwise
            f = f * (1 + noise.level*r);
        end
    end
end

if isfield(options, 'dnoise')
    dnoise = options.dnoise;
    if isstruct(dnoise) && isfield(dnoise, 'level') && dnoise.level > 0
        phi0 = 0.6*cos(1e8*norm(x,9)) + 0.3*sin(100*norm(x,1))*cos(100*norm(x,Inf)) + 0.1*cos(norm(x));
        noisimul = phi0*(4*phi0^2-3);
        switch lower(dnoise.type)
        case {'absolute', 'additive', 'add', 'a', '+'}
            f = f + dnoise.level*noisimul;
        otherwise
            f = f * (1 + dnoise.level*noisimul);
        end
    end
end

if isfield(options, 'single') && isscalar(options.single) && islogical(options.single) && options.single
    f = double(single(f));
end

if (isfield(options, 'signif'))
    sig = min(max(1, options.signif), 16);
    sf = eval(mat2str(f, sig));
    r = sin(sin(sig) + sin(1e8*f) + sum(abs(sin(1e8*x))) + sin(length(x)));
    f = sf + (f-sf)*(r+1);   % This makes the truncation more "irregular".
end

return


function f = evalfun(fun, x, options, ir)
if isfield(options, 'single') && isscalar(options.single) && islogical(options.single) && options.single
    f = fun(double(single(x)));
else
    f = fun(x);
end
options.ir = ir;
f = evalf(f, x, options);
return


function [cineq, ceq] = evalcon(con, x, options, ir)
if isfield(options, 'single') && isscalar(options.single) && islogical(options.single) && options.single
    [cineq, ceq] = con(double(single(x)));
else
    [cineq, ceq] = con(x);
end
options.ir = ir;
afun = @(f) evalf(f, x, options);
cineq = arrayfun(afun, cineq);
ceq = arrayfun(afun, ceq);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cstrv = get_cstrv(x, Aineq, bineq, Aeq, beq, lb, ub, nonlcon)
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
if ~isempty(nonlcon)
    [nlcineq, nlceq] = nonlcon(x);
else
    nlcineq = [];
    nlceq = [];
end
cstrv = max([0; rineq; abs(req); lb-x; x-ub; nlcineq; abs(nlceq)], [], 'includenan');
% max(X, [], 'includenan') returns NaN if X contains NaN, and maximum of X otherwise
return
