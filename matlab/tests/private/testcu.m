function [mrec, output] = testcu(solvers, options)

% Default options
rhobeg = 1;
rhoend = 1e-6;
maxfun_dim = 500;
maxfun = 50000;
maxit = 1000;
ftarget = -inf;
randomizex0 = 0;
noiselevel = 0;
dnoiselevel = 0;
nr = 5;
ctol = 1e-10;
cpenalty = 1e10;
type = 'ubln'; % The types of problems to test
mindim = 1; % The minimal dimension of problems to test
if strcmpi(solvers{1}, 'cobyla') || strcmpi(solvers{2}, 'cobyla')
    maxdim = 20; % The maximal dimension of problems to test
else
    maxdim = 50; % The maximal dimension of problems to test
end
mincon = 0; % The minimal number of constraints of problems to test
maxcon = min(5000, 100*maxdim); % The maximal number of constraints of problems to test
sequential = false;
thorough_test = 0;

% Set options
options = setopt(options, rhobeg, rhoend, maxfun_dim, maxfun, maxit, ftarget, randomizex0, ...
    noiselevel, dnoiselevel, nr, ctol, cpenalty, type, mindim, maxdim, mincon, maxcon, sequential, thorough_test);

% Select the problems to test.
requirements = struct();
requirements.mindim = options.mindim;
requirements.maxdim = options.maxdim;
requirements.mincon = options.mincon;
requirements.maxcon = options.maxcon;
requirements.type = options.type;
if strcmpi(solvers{1}, 'cobyla') || strcmpi(solvers{2}, 'cobyla')
    requirements.blacklist = {};
    requirements.blacklist = [requirements.blacklist, {'CHEBYQADNE','HAIFAM','HIMMELBI','HYDCAR20','LUKSAN12','LUKSAN13','MSS1','SPANHYD','VANDERM1','VANDERM2','VANDERM3', 'TAX13322', 'TAXR13322'}]; % Takes more than 2 min to solve
    requirements.blacklist = [requirements.blacklist, {'DMN15102', 'DMN15103', 'DMN15332', 'DMN15333', 'DMN37142', 'DMN37143'}]; % Time-consuming
end
plist = secup(requirements);

np = length(plist);
ns = length(solvers);
nr = options.randrun;
maxfun = options.maxfun;
sequential = options.sequential;

% These arrays will record the function values and constraint values during the tests.
frec = NaN(np, ns, nr, maxfun);
crec = NaN(np, ns, nr, maxfun);
pdim = NaN(np, 1);  % Data profile needs the dimension of the problem.

if sequential
    for ip = 1 : np
        orig_warning_state = warnoff(solvers);

        pname = plist{ip};

        fprintf('\n%3d. \t%s:\n', ip, upper(pname));

        prob = macup(pname);
        pdim(ip) = length(prob.x0);

        for ir = 1 : nr
            for is = 1 : ns
                [frec(ip, is, ir, :), crec(ip, is, ir, :)] = testsolv(solvers{is}, prob, options);
            end
        end

        warning(orig_warning_state); % Restore the behavior of displaying warnings
    end
else
    parfor ip = 1 : np
        orig_warning_state = warnoff(solvers);

        pname = plist{ip};

        fprintf('\n%3d. \t%s:\n', ip, upper(pname));

        prob = macup(pname);
        pdim(ip) = length(prob.x0);

        for ir = 1 : nr
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

output = struct();
output.plist = plist;
output.pdim = pdim;

return


function [fval_history, cv_history] = testsolv(solver, prob, options)

prob.options = setsolvopt(solver, length(prob.x0), options); % Set the options for the solver

if ischstr(solver)
    if any(regexp(solver, '_classical'))
        prob.options.classical = true;
        solver = regexprep(solver, '_classical', '');
    end
    solver = str2func(solver);
end

maxfun = options.maxfun;
fval_history = NaN(maxfun, 1);
cv_history = NaN(maxfun, 1);

[~, ~, ~, output] = solver(prob);

nf = min(maxfun, output.funcCount); % The solvers may not respect maxfun.

if (nf >= 1)
    fval_history(1:nf) = output.fhist(1:nf);
    if isfield(output, 'chist')
        cv_history(1:nf) = max(0, output.chist(1:nf));
    else
        cv_history(1:nf) = zeros(nf, 1);
    end
    fval_history(nf+1:maxfun) = fval_history(nf);
    cv_history(nf+1:maxfun) = cv_history(nf);
else
    % Sometimes pdfo may return nf = 0, e.g., when it detects infeasiblity.
    fval_history = prob.f0;
    cv_history = prob.constrv0;
end

return


function options = setopt(options, rhobeg, rhoend, maxfun_dim, maxfun, maxit, ftarget, randomizex0, ...
        noiselevel, dnoiselevel, nr, ctol, cpenalty, type, mindim, maxdim, mincon, maxcon, sequential, thorough_test) % Set options

if (~isfield(options, 'rhoend'))
    options.rhoend = rhoend;
end
if (~isfield(options, 'rhobeg'))
    options.rhobeg = rhobeg;
end
if (~isfield(options, 'maxfun_dim'))
    options.maxfun_dim = maxfun_dim;
end
if (~isfield(options, 'maxfun'))
    options.maxfun = maxfun;
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
if (~isfield(options, 'randomizex0'))
    options.randomizex0 = randomizex0;
end
if (~isfield(options, 'noiselevel'))
    options.noiselevel = noiselevel;
end
if (~isfield(options, 'dnoiselevel'))
    options.dnoiselevel = dnoiselevel;
end
if (~isfield(options, 'randrun'))
    options.randrun = nr;
end
if (options.randomizex0 == 0 && options.noiselevel == 0)
    options.randrun = 1;
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
if (~isfield(options, 'sequential'))
    options.sequential = sequential;
end
if (~isfield(options, 'thorough_test'))
    options.thorough_test = thorough_test;
end

return



function solv_options = setsolvopt(solv, n, options)

solv_options = struct();
solv_options.rhobeg = options.rhobeg;
solv_options.rhoend = options.rhoend;
solv_options.maxfun = min(options.maxfun_dim*n, options.maxfun);
solv_options.ftarget = options.ftarget;
solv_options.output_xhist = false;
solv_options.output_nlchist = false;
solv_options.iprint = 0;
solv_options.quiet = true;
solv_options.debug = false;
solv_options.chkfunval = false;
solv_options.classical = false;
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
