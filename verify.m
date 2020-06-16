function success = verify(solvers, options)

success = true;

if nargin < 1
    fprintf('\nSolvers must be specified.\n');
    return
end

if length(solvers) ~= 2 
    fprintf('\nThere should be two solvers.\n')
    return
end

if nargin == 1
    options = struct();
end

if isfield(options, 'prec')
    prec = options.prec;
else
    prec = 0;
end
if isfield(options, 'nr')
    nr = options.nr;
else
    nr = 10;
end

test_options = struct();
test_options.debug = true;
test_options.chkfunval = true;

requirements = struct();
if (isfield(options, 'mindim'))
    requirements.mindim = options.mindim;
else
    requirements.mindim = 1;
end
if (isfield(options, 'maxdim'))
    requirements.maxdim = options.maxdim;
else
    requirements.maxdim = 20;
end
if (isfield(options, 'mincon'))
    requirements.mincon = options.mincon;
else
    requirements.mincon = 0;
end
if (isfield(options, 'maxcon'))
    requirements.maxcon = options.maxcon;
else
    requirements.maxcon = 100;
end
if (isfield(options, 'type'))
    requirements.type = options.type;
else
    requirements.type = 'u';
%    requirements.type = 'ubln';
end

% Supress the following warning
orig_warning_state = warning;
cellfun(@(solver) warning('off', [solver, ':Debug']), solvers);
cellfun(@(solver) warning('off', [solver, ':ChkFunval']), solvers);
cellfun(@(solver) warning('off', [solver, ':Classical']), solvers);
cellfun(@(solver) warning('off', [solver, ':ReviseX0']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownProbField']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownOption']), solvers);

plist = secup(requirements);

fprintf('\n')
for ip = 1 : length(plist)
    pname = plist{ip};
    fprintf('%3d. \t%16s:\t', ip, pname);
    prob = macup(pname);
    x0 = prob.x0;
    n = length(x0);
    for ir = 0 : nr + 4 
        % Some randomization
        rng(ceil(1e5*abs(sin(1e10*(ir+nr)))));
        prob.x0 = x0 + 0.5*randn(size(x0));
        test_options = struct();
        test_options.rhobeg = 1 + 0.5*(2*rand-1);
        test_options.rhoend = 1e-3*(1 + 0.5*(2*rand-1));
        test_options.npt = max(min(ceil(10*rand*n + 2), (n+2)*(n+1)/2), n+2);
        test_options.maxfun = max(ceil(20*n*(1+rand)), n+3);
        test_options.ftarget = -inf;
        if ir == 0
            test_options.npt = (n+2)*(n+1)/2;
        end
        if ir == 1 
            test_options.npt = n + 2;
        end
        if ir == 2 
            test_options.maxfun = test_options.npt + 1;
        end
        if ir == 3
            test_options.rhoend = test_options.rhobeg;
        end 
        if ir == 4
            test_options.ftarget = inf;
        end 
        prob.options = test_options;
        [x1, fx1, exitflag1, output1] = feval(solvers{1}, prob);
        if output1.funcCount == test_options.maxfun && exitflag1 == 0
            exitflag1 = 3;
            display('exitflag1 changed from 0 to 3.')
        end
        display('++++++++')
        [x2, fx2, exitflag2, output2] = feval(solvers{2}, prob);
        if ~iseq(x1, fx1, exitflag1, output1, x2, fx2, exitflag2, output2, prec)
            fprintf('The solvers produce different results on %s at the %dth run.\n', pname, ir);
            success = false;
            keyboard
        end
    end
    decup(pname);
    if success
        fprintf('Success\n');
    else
        fprintf('FAIL!\n')
    end
end

warning(orig_warning_state); % Restore the behavior of displaying warnings

return


function eq = iseq(x, f, exitflag, output, xx, ff, ee, oo, prec) 
eq = true;

if ~isempty(setdiff(fieldnames(output), fieldnames(oo))) || ~isempty(setdiff(fieldnames(oo), fieldnames(output)))
    eq = false;
end

if ~isfield(output,'constrviolation')
    output.constrviolation = 0;
end
if ~isfield(oo,'constrviolation')
    oo.constrviolation = 0;
end

if ~isfield(output, 'chist')
    output.chist = zeros(output.funcCount);
end
if ~isfield(oo, 'chist')
    oo.chist = zeros(oo.funcCount);
end

if (norm(xx-x)/(1+norm(x)) > prec || abs(ff-f)/(1+abs(f)) > prec || abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation)) > prec)
    eq = false;
end

if length(output.fhist) ~= length(oo.fhist) || norm(output.fhist-oo.fhist)/(1+norm(output.fhist)) > prec
    eq = false;
end

if length(output.chist) ~= length(oo.chist) || norm(output.chist-oo.chist)/(1+norm(output.chist)) > prec
    eq = false;
end

if (prec == 0 && (exitflag ~= ee|| oo.funcCount ~= output.funcCount))
    eq = false;
end

diff = max([abs(ff-f)/(1+abs(f)), norm(xx-x)/(1+norm(x)), abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation))]); 

return
