function success = verify(solvers, options)

% Last Modified: Mon 24 May 2021 04:17:01 PM HKT

% As an example:
% options=[]; options.maxdi=20; options.nr=20; verify({'newuoan', 'newuoa'}, options)

% verifies newuoan against newuoa on problems of at most 20 variables, 20 random runs for each
% problem.
% NOTE that newuoa has to be the version in OPDFO, which has been modified (slightly) to behave
% the same as newuoan.

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
    nr = 20;
end
if isfield(options, 'ir')
    % ir is the index of the random experiment to be conducted. If it is negative, then experiments
    % 1, ..., nr, ..., nr + 20 will be conducted. nr + 20 is because there are a few fixed
    % experiments that will always be run.
    ir = options.ir;
else
    ir = -1;
end

if isfield(options, 'minip')
    minip=options.minip;
else
    minip = 1;
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
    requirements.maxdim = 40;
end
if (isfield(options, 'mincon'))
    requirements.mincon = options.mincon;
else
    requirements.mincon = 0;
end
if (isfield(options, 'maxcon'))
    requirements.maxcon = options.maxcon;
else
    requirements.maxcon = min(5000, 100*requirements.maxdim);
end
if (isfield(options, 'type'))
    requirements.type = options.type;
else
    requirements.type = 'ubln';
end

% Suppress the following warning
orig_warning_state = warning;
cellfun(@(solver) warning('off', [solver, ':Debug']), solvers);
cellfun(@(solver) warning('off', [solver, ':ChkFunval']), solvers);
cellfun(@(solver) warning('off', [solver, ':ReviseX0']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownProbField']), solvers);
cellfun(@(solver) warning('off', [solver, ':UnknownOption']), solvers);
cellfun(@(solver) warning('off', [solver, ':InvalidMaxfun']), solvers);
cellfun(@(solver) warning('off', [solver, ':ExtremeBarrier']), solvers);
cellfun(@(solver) warning('off', [solver, ':IprintContradictFortran']), solvers);
cellfun(@(solver) warning('off', [solver, ':InvalidMaxhist']), solvers);
cellfun(@(solver) warning('off', [solver, ':InvalidNpt']), solvers);

if isempty(requirements.list)
    plist = secup(requirements);
else
    plist = requirements.list; % Use the list provided by the user
end

if ir < 0
    minir = 0;
    maxir = nr + 20;
else
    minir = ir;
    maxir = ir;
end

blacklist = [];
%blacklist={'gauss2', 'gauss3','HS25NE', 'cubene'};  % Memory error
blacklist=[blacklist, 'qcnew'];  % f(x) ~= fx due to bad condition
blacklist=[blacklist, 'BLEACHNG']; % takes too much time

fprintf('\n')
for ip = minip : length(plist)
    pname = upper(plist{ip});
    if ismember(lower(pname), lower(blacklist))
        fprintf ('\n!!!!!! %s is skipped!!!!!!\n', pname)
        continue
    end
    fprintf('\n%3d. \t%16s:\t', ip, pname);
    prob = macup(pname);
    objective = prob.objective;
    nonlcon = prob.nonlcon;
    x0 = prob.x0;
    n = length(x0);
    fprintf('\n');
    for ir = minir : maxir
        fprintf('\n%s Run No. %3d: \t', pname, ir);
        % Some randomization
        rng(ceil(1e6*abs(cos(1e6*sin(1e6*(sum(double(pname))*n*ir))))));
        prob.x0 = x0 + 0.5*randn(size(x0));
        test_options = struct();
        test_options.debug = true;
        test_options.chkfunval = true;
        test_options.rhobeg = 1 + 0.5*(2*rand-1);
        test_options.rhoend = 1e-3*(1 + 0.5*(2*rand-1));
        test_options.npt = max(min(ceil(10*rand*n + 2), (n+2)*(n+1)/2), n+2);
        test_options.maxfun = max(ceil(20*n*(1+rand)), n+3);
        test_options.ftarget = -inf;
        %test_options.classical = (randn < -1.2);
        test_options.classical = 0;
        test_options.fortran = (rand > 0.5);
        test_options.output_xhist = (rand > 0.5);
        test_options.maxhist = ceil(randn*1.5*test_options.maxfun);
        %test_options.iprint = floor(3*rand);
        %test_options.quiet = (rand > 0.5);
        test_options.iprint = 0;
        test_options.quiet = true;
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
            test_options.maxfun = 1;
        end
        if ir == 5
            test_options.maxfun = ceil(n/2);
        end
        if ir == 6
            test_options.ftarget = inf;
        end
        if ir == 7
            test_options.rhoend = test_options.rhobeg;
        end
        if 0 <= ir && ir <= 20
            test_options.chkfunval = false;  % The checking would fail due to noise.
            prob.objective = @(x) noisyfeval(objective, x);
            if ~isempty(nonlcon)
                prob.nonlcon = @(x) noisyceval(nonlcon, x);
            end
        else
            prob.objective  = objective;
            prob.nonlcon = nonlcon;
        end

        prob.options = test_options;


        %tic;
        [x1, fx1, exitflag1, output1] = feval(solvers{1}, prob);
        %T = toc; fprintf('\nRunning time for %s:\t %f\n', solvers{1}, T);

        %tic;
        [x2, fx2, exitflag2, output2] = feval(solvers{2}, prob);
        %T = toc; fprintf('\nRunning time for %s:\t %f\n', solvers{2}, T);


        if output1.funcCount == test_options.maxfun && (exitflag1 == 0 || exitflag1 == 2) && exitflag2 == 3
            exitflag1 = 3;
            %display('exitflag1 changed to 3.')
        end
        if output2.funcCount == test_options.maxfun && (exitflag2 == 0 || exitflag2 == 2) && exitflag1 == 3
            exitflag2 = 3;
            %display('exitflag2 changed to 3.')
        end
        if fx1 <= test_options.ftarget
            exitflag1 = 1;
            fprintf('\nexitflag1 changed to 1.\n')
        end
        if fx2 <= test_options.ftarget
            exitflag2 = 1;
            fprintf('\nexitflag2 changed to 1.\n')
        end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Special Treatments%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       minfhist = min(length(output1.fhist), length(output2.fhist));
       % NEWUOA
       if (strcmpi(solvers{1}, 'newuoa') || strcmpi(solvers{2}, 'newuoan')) && exitflag1 == 2 && exitflag2 ~=2 && output1.funcCount <= output2.funcCount && all(output2.fhist(end-minfhist+1:end-(output2.funcCount-output1.funcCount)) == output1.fhist(end-minfhist+(output2.funcCount-output1.funcCount)+1:end)) && fx2 <= fx1
           x2 = x1;
           fx2 = fx1;
           exitflag2 = exitflag1;
           output2.fhist = output1.fhist;
           output2.funcCount = output1.funcCount;
           fprintf('\nThe original solver exits due to failure of the TR subproblem solver.\n');
       end
       if (strcmpi(solvers{1}, 'newuoan') || strcmpi(solvers{2}, 'newuoa')) && exitflag2 == 2 && exitflag1 ~=2 &&  output2.funcCount <= output1.funcCount && all(output1.fhist(end-minfhist+1:end-(output1.funcCount-output2.funcCount)) == output2.fhist(end-minfhist+(output1.funcCount-output2.funcCount)+1:end)) && fx1 <= fx2
           x1 = x2;
           fx1 = fx2;
           exitflag1 = exitflag2;
           output1.fhist = output2.fhist;
           output1.funcCount = output2.funcCount;
           fprintf('\nThe original solver exits due to failure of the TR subproblem solver.\n');
       end
       if (strcmpi(solvers{1}, 'newuoa') || strcmpi(solvers{2}, 'newuoa')) && fx1 == fx2 && output1.funcCount == output2.funcCount && all(output1.fhist(end-minfhist+1:end) == output2.fhist(end-minfhist+1:end)) && norm(x1 - x2) > 0
           x1 = x2;
           fprintf('\nx1 changed to x2\n');
       end

       % COBYLA
        if ((strcmpi(solvers{1}, 'cobyla') || strcmpi(solvers{2}, 'cobyla'))  && fx1 == fx2 && (~isfield(output1,'constrviolation') && ~isfield(output2, 'constrviolation') || isfield(output1, 'constrviolation') && output1.constrviolation == output2.constrviolation) && norm(x1-x2)>0)
            x1 = x2;
            fprintf('\nx1 changed to x2.\n');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~iseq(x1, fx1, exitflag1, output1, x2, fx2, exitflag2, output2, prec)
            keyboard
            error('The solvers produce different results on %s at the %dth run.\n', pname, ir);
        end
    end
    decup(pname);
end

if ~isempty(intersect(plist, blacklist))
    fprintf ('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    fprintf ('\n!!!!!! The following problems are skipped!!!!!!\n');
    intersect(plist, blacklist)
    fprintf ('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
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

function f = noisy(f, x, noise_level)
    if nargin < 3
        noise_level = 5e-1;
    end
    r = sin((abs(f)+1)*1e6*max(abs(x)+1)*cos(1e6*sum(abs(x)+1)*sin(1e6*(norm(x)+1))));
    if (r > 0.75)
        r = sign(f)*inf;
    elseif (r > 0.5)
        r = NaN;
    elseif (r < - 0.75)
        r = -sign(f)*inf;
    end
    f = f*(1+noise_level*r);
return

function f = noisyfeval(func, x, noise_level)
    if nargin < 3
        noise_level = 5e-1;
    end
    f = feval(func, x);
    f = noisy(f, x, noise_level);
return

function [cineq, ceq] = noisyceval(con, x, noise_level)
    if nargin < 3
        noise_level = 5e-1;
    end
    [cineq, ceq] = feval(con, x);
    for i = 1 : length(cineq)
        cineq(i) = noisy(cineq(i), x, noise_level);
    end
    for i = 1 : length(ceq)
        ceq(i) = noisy(ceq(i), x, noise_level);
    end
return

function [x,fx, exitflag, output] = newuoan1(varargin)
    [x,fx, exitflag, output] = newuoan(varargin, 1);
return
