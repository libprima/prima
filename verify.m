function success = verify(solvers, options)

% As an example:
% options=[]; options.maxdi=20; options.nr=20; verify({'newuoan', 'newuoa'}, options)

% verifies newuoan against newuoa on problems of at most 20 variables, 20 random runs for each
% problem.
% NOTE that newuoa has to be the version in OPDFO, which has been modified (slightly) to behave
% the same as newuoan.
%
% Last Modified: Monday, May 24, 2021 PM02:19:45 

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
if isfield(options, 'ir')
    % ir is the index of the random experiment to be conducted. If it is negative, then experiments
    % 0, 1, ..., nr, ..., nr + 5 will be conducted. nr + 5 is because there are a few fixed
    % experiments that will always be run.
    ir = options.ir;
else
    ir = -1;
end

requirements = struct();
if isfield(options, 'list')
    requirements.list = options.list;  % Only test problems in this list
else
    requirements.list = {};  % No restriction
end
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
cellfun(@(solver) warning('off', [solver, ':ReviseX0']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownProbField']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownOption']), solvers);
cellfun(@(solver) warning('off', [solver, ':InvalidMaxfun']), solvers);
cellfun(@(solver) warning('off', [solver, ':ExtremeBarrier']), solvers);

if isempty(requirements.list)
    plist = secup(requirements);
else
    plist = requirements.list; % Use the list provided by the user
end

if ir < 0
    minir = 0;
    maxir = nr + 5;
else
    minir = ir;
    maxir = ir;
end

fprintf('\n')
for ip = 1 : length(plist)
    pname = plist{ip};
    fprintf('%3d. \t%16s:\t', ip, pname);
%    fprintf('%3d. \t%16s:\t\n', ip, pname);
    prob = macup(pname);
    x0 = prob.x0;
    n = length(x0);
    if n > 30
        fprintf('\n');
    end
    for ir = minir : maxir
        if n > 30
            fprintf('Run No. %3d: \t', ir);
        end
        % Some randomization
        rng(ceil(1e6*abs(sin(1e6*(sum(double(pname))*n*ip*ir*nr*requirements.mindim*requirements.maxdim))))); %!!! TO MUCH RANDOMNESS!!! DIFFICULT TO REPRODUCE!!!
        %rng(ceil(1e6*abs(sin(1e6*(sum(double(pname))*n*ip*ir)))));
        prob.x0 = x0 + 0.5*randn(size(x0));
        test_options = struct();
        test_options.debug = true;
        test_options.chkfunval = true;
        test_options.rhobeg = 1 + 0.5*(2*rand-1);
        test_options.rhoend = 1e-3*(1 + 0.5*(2*rand-1));
        test_options.npt = max(min(ceil(10*rand*n + 2), (n+2)*(n+1)/2), n+2);
        test_options.maxfun = max(ceil(20*n*(1+rand)), n+3);
        test_options.ftarget = -inf;
        test_options.classical = (randn < -1.2);
        test_options.fortran = (rand > 0.5);
        test_options.output_xhist = (rand > 0.5);
        test_options.maxhist = ceil(randn*1.5*test_options.maxfun);
        test_options.iprint = floor(3*rand);
        test_options.quiet = (rand > 0.5);
        if mod(ir, 50) == 0 && exist('NEWUOA_output.txt', 'file')
            delete('NEWUOA_output.txt');
        end
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
            test_options.maxfun = 1000*n;
        end
        if ir == 4
            test_options.ftarget = inf;
        end
        if ir == 5
            test_options.rhoend = test_options.rhobeg;
        end
        prob.options = test_options;

        tic;
        [x1, fx1, exitflag1, output1] = feval(solvers{1}, prob);
        T = toc;
        fprintf('\nRunning time for %s:\t %f\n', solvers{1}, T);

        tic;
        [x2, fx2, exitflag2, output2] = feval(solvers{2}, prob);
        T = toc;
        fprintf('\nRunning time for %s:\t %f\n', solvers{2}, T);


        if output1.funcCount == test_options.maxfun && (exitflag1 == 0 || exitflag1 == 2) && exitflag2 == 3
            exitflag1 = 3;
            %display('exitflag1 changed to 3.')
        end
        if output2.funcCount == test_options.maxfun && (exitflag2 == 0 || exitflag2 == 2) && exitflag1 == 3
            exitflag2 = 3;
            %display('exitflag2 changed to 3.')
        end
        if iseq(x1, fx1, exitflag1, output1, x2, fx2, exitflag2, output2, prec)
            if n > 30
                fprintf('Succeed\n');
            end
        else
            fprintf('The solvers produce different results on %s at the %dth run.\n', pname, ir);
            success = false;
            keyboard
        end
    end
    decup(pname);
    if success
        fprintf('Succeed\n');
    else
        fprintf('FAIL!\n')
    end
end

warning(orig_warning_state); % Restore the behavior of displaying warnings

return


function eq = iseq(x, f, exitflag, output, xx, ff, ee, oo, prec)
eq = true;

if ~isempty(setdiff(fieldnames(output), [fieldnames(oo); 'fhist'; 'xhist'])) || ~isempty(setdiff(fieldnames(oo), [fieldnames(output); 'fhist'; 'xhist']))
    eq = false;
end

if ~isfield(output,'constrviolation')
    output.constrviolation = 0;
end
if ~isfield(oo,'constrviolation')
    oo.constrviolation = 0;
end

if ~isfield(output, 'chist')
    output.chist = zeros(output.funcCount, 1);
end
if ~isfield(oo, 'chist')
    oo.chist = zeros(oo.funcCount, 1);
end

if (norm(xx-x)/(1+norm(x)) > prec || abs(ff-f)/(1+abs(f)) > prec || abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation)) > prec)
    eq = false;
end

if isfield(output, 'fhist')
    output.fhist = output.fhist(:);
else
    output.fhist = [];
end
if isfield(oo, 'fhist')
    oo.fhist = oo.fhist(:);
else
    oo.fhist = [];
end
nhist = min(length(output.fhist), length(oo.fhist));
output.fhist = output.fhist(end - nhist + 1: end);
oo.fhist = oo.fhist(end - nhist + 1: end);

if norm(output.fhist(1:nhist)-oo.fhist(1:nhist))/(1+norm(output.fhist(1:nhist))) > prec
    eq = false;
end

if norm(output.chist(1:nhist)-oo.chist(1:nhist))/(1+norm(output.chist(1:nhist))) > prec
    eq = false;
end

if (prec == 0 && (exitflag ~= ee|| oo.funcCount ~= output.funcCount))
    eq = false;
end

diff = max([abs(ff-f)/(1+abs(f)), norm(xx-x)/(1+norm(x)), abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation))]);

return
