function isequiv(solvers, options)
%This function verifies that two solvers produce the same result on CUTEst problems.
%
% As an example:
% options=[]; options.maxdi=20; options.nr=20; isequiv({'newuoa', 'newuoa_norma'}, options)
%
% verifies newuoa against newuoa_norma on problems of at most 20 variables, 20 random runs for each problem.
%
% Coded by Zaikun ZHANG (www.zhangzk.net).
%
% Started: July 2020
%
% Last Modified: Monday, April 07, 2023 PM07:04:00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 1
    fprintf('\nSolvers must be specified.\n');
    return
end

if length(solvers) ~= 2
    fprintf('\nThere should be two solvers.\n')
    return
end

solvers = lower(solvers);

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
% ir is the index of the random experiment to be conducted. If it is negative, then experiments
% 1, ..., nr, ..., nr + 20 will be conducted. nr + 20 is because there are fixed experiments
% that will always be run.
if isfield(options, 'ir')
    ir = options.ir;
else
    ir = -1;
end

if ir < 0
    if isfield(options, 'minir')
        minir = options.minir;
    else
        minir = 0;
    end
    if isfield(options, 'maxir')
        maxir = options.maxir;
    else
        maxir = nr + 20;
    end
else
    minir = ir;
    maxir = ir;
end

if isfield(options, 'minip')
    minip=options.minip;
else
    minip = 1;
end

if isfield(options, 'maxip')
    maxip=options.maxip;
else
    maxip = 2^32 - 1;
end

requirements = struct();
if isfield(options, 'list')
    requirements.list = options.list;  % Only test problems in this list
else
    requirements.list = {};  % No restriction
end
if isfield(options, 'blacklist')
    requirements.blacklist = options.blacklist;
else
    requirements.blacklist = {};
end
if (isfield(options, 'mindim'))
    requirements.mindim = options.mindim;
else
    requirements.mindim = 1;
end
if (isfield(options, 'maxdim'))
    requirements.maxdim = options.maxdim;
else
    if any(startsWith(lower(solvers), 'cobyla'))
        requirements.maxdim = 20;
    else
        requirements.maxdim = 50;
    end
end
if (isfield(options, 'mincon'))
    requirements.mincon = options.mincon;
else
    requirements.mincon = 0;
end
if (isfield(options, 'maxcon'))
    requirements.maxcon = options.maxcon;
else
    requirements.maxcon = min(20000, 100*requirements.maxdim);
end
if (isfield(options, 'type'))
    requirements.type = options.type;
else
    requirements.type = 'ubln';
end

if ~isempty(requirements.list)
    plist = requirements.list; % Use the list provided by the user, neglecting all other requirements
    if (ischarstr(plist))  % In case plist is indeed the name of a problem
        plist = {plist};
    end
else
    requirements.blacklist = [requirements.blacklist, black_list(solvers{1}), black_list(solvers{2})];
    plist = secup(requirements);
end

np = length(plist);

single_test = (np <= 1);
if isfield(options, 'sequential')
    sequential = options.sequential;
else
    sequential = single_test;
end

maxip = min(np, maxip);

% Directories for recording the starting/ending of problems (tic/toc are unavailable in parfor).
stamp = strjoin(solvers, '_');
prob_start_dir = strtrim(fullfile(options.test_dir, [stamp, '_start']));
prob_start_time_dir = strtrim(fullfile(options.test_dir, [stamp, '_start_time']));
prob_start_runs_dir = strtrim(fullfile(options.test_dir, [stamp, '_start_runs']));
prob_end_dir = strtrim(fullfile(options.test_dir, [stamp, '_end']));
prob_end_time_dir = strtrim(fullfile(options.test_dir, [stamp, '_end_time']));
prob_end_runs_dir = strtrim(fullfile(options.test_dir, [stamp, '_end_runs']));
system(['rm -rf ', prob_start_dir, '; ', 'mkdir -p ', prob_start_dir]);
system(['rm -rf ', prob_start_time_dir, '; ', 'mkdir -p ', prob_start_time_dir]);
system(['rm -rf ', prob_start_runs_dir, '; ', 'mkdir -p ', prob_start_runs_dir]);
system(['rm -rf ', prob_end_dir, '; ', 'mkdir -p ', prob_end_dir]);
system(['rm -rf ', prob_end_time_dir, '; ', 'mkdir -p ', prob_end_time_dir]);
system(['rm -rf ', prob_end_runs_dir, '; ', 'mkdir -p ', prob_end_runs_dir]);
disp(['prob_start_dir = ', prob_start_dir]);
disp(['prob_start_time_dir = ', prob_start_time_dir]);
disp(['prob_start_runs_dir = ', prob_start_runs_dir]);
disp(['prob_end_dir = ', prob_end_dir]);
disp(['prob_end_time_dir = ', prob_end_time_dir]);
disp(['prob_end_runs_dir = ', prob_end_runs_dir]);

if sequential
    for ip = minip : maxip

        % Turn off unwanted warnings
        orig_warning_state = warnoff(solvers);

        pname = upper(plist{ip});
        [~, time] = system('date +%y%m%d_%H%M%S');
        system(['touch ', fullfile(prob_start_dir, pname)]);
        system(['touch ', fullfile(prob_start_time_dir, [pname, '.', strtrim(time)])]);
        system(['touch ', fullfile(prob_start_runs_dir, pname)]);
        fprintf('\n%3d. \t%s starts at %s\n', ip, pname, char(datetime()));

        prob = macup(pname);

        for ir = minir : maxir
            % The following line compares the solvers on `prob`; ir is needed for the random seed, and
            % `prec` is the precision of the comparison (should be 0). The function will raise an error
            % if the solvers behave differently.
            system(['touch ', fullfile(prob_start_runs_dir, [pname, '.', num2str(ir, '%02d')])]);
            fprintf('\n     \t%s Run No. %3d starts at %s\n', pname, ir, char(datetime()));
            compare(solvers, prob, ir, prec, single_test, options);
            system(['touch ', fullfile(prob_end_runs_dir, [pname, '.', num2str(ir, '%02d')])]);
            fprintf('\n     \t%s Run No. %3d ends at %s\n', pname, ir, char(datetime()));
        end

        decup(prob);

        [~, time] = system('date +%y%m%d_%H%M%S');
        system(['touch ', fullfile(prob_end_dir, pname)]);
        system(['touch ', fullfile(prob_end_time_dir, [pname, '.', strtrim(time)])]);
        fprintf('\n%3d. \t%s ends at %s\n', ip, pname, char(datetime()));

        % Restore the behavior of displaying warnings
        warning(orig_warning_state);
    end
else
    parfor ip = minip : maxip

        % Turn off unwanted warnings
        orig_warning_state = warnoff(solvers);

        pname = upper(plist{ip});
        [~, time] = system('date +%y%m%d_%H%M%S');
        system(['touch ', fullfile(prob_start_dir, pname)]);
        system(['touch ', fullfile(prob_start_time_dir, [pname, '.', strtrim(time)])]);
        fprintf('\n%3d. \t%s starts at %s\n', ip, pname, char(datetime()));

        prob = macup(pname);

        for ir = minir : maxir
            % The following line compares the solvers on `prob`; ir is needed for the random seed, and
            % `prec` is the precision of the comparison (should be 0). The function will raise an error
            % if the solvers behave differently.
            system(['touch ', fullfile(prob_start_runs_dir, [pname, '.', num2str(ir, '%02d')])]);
            fprintf('\n     \t%s Run No. %3d starts at %s\n', pname, ir, char(datetime()));
            compare(solvers, prob, ir, prec, single_test, options);
            system(['touch ', fullfile(prob_end_runs_dir, [pname, '.', num2str(ir, '%02d')])]);
            fprintf('\n     \t%s Run No. %3d ends at %s\n', pname, ir, char(datetime()));
        end

        decup(prob);

        [~, time] = system('date +%y%m%d_%H%M%S');
        system(['touch ', fullfile(prob_end_dir, pname)]);
        system(['touch ', fullfile(prob_end_time_dir, [pname, '.', strtrim(time)])]);
        fprintf('\n%3d. \t%s ends at %s\n', ip, pname, char(datetime()));

        % Restore the behavior of displaying warnings
        warning(orig_warning_state);
    end
end

disp(['prob_start_dir = ', prob_start_dir]);
disp(['prob_start_time_dir = ', prob_start_time_dir]);
disp(['prob_start_runs_dir = ', prob_start_runs_dir]);
disp(['prob_end_dir = ', prob_end_dir]);
disp(['prob_end_time_dir = ', prob_end_time_dir]);
disp(['prob_end_runs_dir = ', prob_end_runs_dir]);

fprintf('\n\nSucceed!\n\n');   % Declare success if we arrive here without an error.

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function equiv = compare(solvers, prob, ir, prec, single_test, options)
pname = prob.name;
objective = prob.objective;
nonlcon = prob.nonlcon;
x0 = prob.x0;
n = length(x0);

% Some randomization
% Set seed using pname, n, and ir. We ALTER THE SEED weekly to test the solvers as much as possible.
% N.B.: The weeknum function considers the week containing January 1 to be the first week of the
% year, and increments the number every SUNDAY.
if isfield(options, 'yw')
    yw = options.yw;
elseif isfield(options, 'seed')
    yw = options.seed;
else
    yw = year_week('Asia/Shanghai');
end
fprintf('\nYW = %d\n', yw);
rseed = max(0, min(2^32 - 1,  sum(pname) + n + ir + yw));  % A random seed defined by the current test and yw
orig_rng_state = rng();  % Save the current random number generator settings
rng(rseed);  % Set the random seed for reproducibility
prob.x0 = x0 + 0.5*randn(size(x0));
test_options = struct();
test_options.rhobeg = 1 + 0.5*(2*rand-1);
test_options.rhoend = 1e-3*(1 + 0.5*(2*rand-1));
test_options.npt = max(min(floor(6*rand*n), (n+2)*(n+1)/2), n+2);
test_options.maxfun = max(ceil(20*n*(1+rand)), n+3);  % For reproducibility, do not remove this even if `options` contains `maxfun`.
if isfield(options, 'maxfun')
    test_options.maxfun = options.maxfun;
end
test_options.ftarget = objective(x0) - 10*abs(randn)*max(1, objective(x0));
%test_options.fortran = (rand > 0.5);
test_options.fortran = true;
test_options.output_xhist = (rand > 0.5);
%test_options.output_xhist = 1;
test_options.output_nlchist = (rand > 0.5);
test_options.maxhist = ceil(randn*1.5*test_options.maxfun);
%test_options.maxhist = test_options.maxfun;
if single_test
    % DO NOT INVOKE ANY RANDOMIZATION WITHIN THIS IF. Otherwise, a single test cannot reproduce the
    % corresponding test in a multiple one.
    test_options.maxhist = test_options.maxfun;
    test_options.output_xhist = true;
    test_options.output_nlchist = true;
end
test_options.maxfilt = ceil(randn*500);
test_options.iprint = floor(3*rand);
test_options.quiet = (rand < 0.9);
% Test all precisions. For unavailable precisions, the double-precision version will be called.
if rand < 0.7  % Prob = 0.6
    test_options.precision = 'double';
elseif rand < 0.9  % Prob = 0.27
    test_options.precision = 'single';
else  % Prob = 0.03
    test_options.precision = 'quadruple';
end

%!------------------------------------------------------------------------------------------------!%
% Test both debugging and non-debugging versions. They may behave differently.
% On 20220302, it is observed that, when the Fortran code is compiled with the '-g' (debugging)
% option, the INTENT(OUT) arguments will keep the values that they get before entering subroutines,
% even though such values should be cleared on entry of the subroutines. This behavior makes it
% impossible to detect the arguments that should be INTENT(INOUT) but mistakenly declared as
% INTENT(OUT). The observation was made on the argument named SNORM in the subroutine TRSTEP of
% LINCOA, and it took a whole day to debug.
test_options.debug = (rand < 0.7);
test_options.chkfunval = test_options.debug;
%!------------------------------------------------------------------------------------------------!%

% Test all variants. If the classical variant is unavailable, the modernized variant will be called.
test_options.classical = ~(isfield(options, 'no_classical') && options.no_classical) && (rand < 0.1);
% Test only double for the classical variant; debugging version is unavailable for the classical variant.
if test_options.classical
    test_options.precision = 'double';
    test_options.debug = false;
    test_options.chkfunval = false;
end

call_by_package = (rand < 0.5);  % Call by the package instead of the solver
call_by_structure = (rand < 0.5);  % Pass the problem by a structure
if mod(ir, 50) == 0 && ~isempty(dir('*_output.txt'))
    delete('*_output.txt');
end
if ir == 1
    test_options.npt = (n+2)*(n+1)/2;
end
if ir == 2
    test_options.npt = n + 2;
end
if ir == 3
    test_options.maxfun = test_options.npt + 1;
end
if ir == 4
    test_options.maxfun = 1000*n;
end
if ir == 5
    test_options.maxfun = 1;
end
if ir == 6
    test_options.maxfun = ceil(n/2);
end
if ir == 7
    test_options.ftarget = inf;
end
if ir == 8
    test_options.rhoend = test_options.rhobeg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ir == 9
    test_options.npt = 2*n;
end
if 10 <= ir && ir <= 12
    test_options.npt = ceil(rand*n^2);
end
if 13 <= ir && ir <= 15
    test_options.npt = floor(2*rand*n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 <= ir && ir <= 20
    % The TOUGH tests
    % We must pass the random seed `rseed` to `tough` to ensure reproducibility.
    test_options.chkfunval = false;  % The checking would fail due to noise.
    prob = tough(prob, rseed);
else
    prob.objective  = objective;
    prob.nonlcon = nonlcon;
end
prob.options = test_options;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN: Call the solvers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N.B.: In some tests, we may invoke this function with solvers{1} == solvers{2}. So do NOT assume
% that one of the solvers is 'SOLVER' and the other is 'SOLVER_norma'.

% Use function handle to avoid `feval`.
solver1 = str2func(solvers{1});
solver2 = str2func(solvers{2});

if endsWith(solvers{1}, '_norma')
    package1 = @prima_norma;
else
    package1 = @prima;
end
if endsWith(solvers{2}, '_norma')
    package2 = @prima_norma;
else
    package2 = @prima;
end

tested_solver_name = regexprep(solvers{1}, '_norma', '');

prob1 = prob;
prob2 = prob;

test_row_x = (rand > 0.5);
test_row_bounds = (rand > 0.5);
if ~endsWith(solvers{1}, '_norma')
    if test_row_x
        prob1.x0 = prob1.x0';
        prob1.objective = @(x) prob1.objective(x');
        if ~isempty(prob1.nonlcon)
            prob1.nonlcon = @(x) prob1.nonlcon(x');
        end
    end
    if test_row_bounds
        prob1.lb = prob1.lb';
        prob1.ub = prob1.ub';
    end
end
if ~endsWith(solvers{2}, '_norma')
    if test_row_x
        prob2.x0 = prob2.x0';
        prob2.objective = @(x) prob2.objective(x');
        if ~isempty(prob2.nonlcon)
            prob2.nonlcon = @(x) prob2.nonlcon(x');
        end
    end
    if test_row_bounds
        prob2.lb = prob2.lb';
        prob2.ub = prob2.ub';
    end
end

test_fixed_x = (rand > 0.5 && ~isempty(prob1.lb) && ~isempty(prob1.ub));
if test_fixed_x
    fixedxl = (rand(n, 1) > 0.6);
    fixedxu = (rand(n, 1) > 0.3 & ~fixedxl);
    prob1.ub(fixedxl) = prob1.lb(fixedxl);
    prob1.lb(fixedxu) = prob1.ub(fixedxu);
    prob2.lb = prob1.lb;
    prob2.ub = prob1.ub;
end

test_infnan_lcon = (rand > 0.5 && ~(isempty(prob1.Aineq) && isempty(prob1.Aeq)));
if test_infnan_lcon
    if ~isempty(prob1.Aineq)
        prob1.Aineq(rand(size(prob1.Aineq)) > 0.9) = Inf;
        prob1.Aineq(rand(size(prob1.Aineq)) > 0.9) = -Inf;
        prob1.Aineq(rand(size(prob1.Aineq)) > 0.9) = NaN;
        prob1.bineq(rand(size(prob1.bineq)) > 0.9) = Inf;
        prob1.bineq(rand(size(prob1.bineq)) > 0.9) = -Inf;
        prob1.bineq(rand(size(prob1.bineq)) > 0.9) = NaN;
    end
    if ~isempty(prob1.Aeq)
        prob1.Aeq(rand(size(prob1.Aeq)) > 0.9) = Inf;
        prob1.Aeq(rand(size(prob1.Aeq)) > 0.9) = -Inf;
        prob1.Aeq(rand(size(prob1.Aeq)) > 0.9) = NaN;
        prob1.beq(rand(size(prob1.beq)) > 0.9) = Inf;
        prob1.beq(rand(size(prob1.beq)) > 0.9) = -Inf;
        prob1.beq(rand(size(prob1.beq)) > 0.9) = NaN;
    end
    prob2.Aineq = prob1.Aineq;
    prob2.bineq = prob1.bineq;
    prob2.Aeq = prob1.Aeq;
    prob2.beq = prob1.beq;
end

test_row_lineq = (rand > 0.5);
test_row_leq = (rand > 0.5);
if ~endsWith(solvers{2}, '_norma')
    if test_row_lineq
        prob1.bineq = prob1.bineq';
    end
    if test_row_leq
        prob1.beq = prob1.beq';
    end
end
if ~endsWith(solvers{2}, '_norma')
    if test_row_lineq
        prob2.bineq = prob2.bineq';
    end
    if test_row_leq
        prob2.beq = prob2.beq';
    end
end

exception = [];
try
    if call_by_package
        if call_by_structure
            prob1.options.solver = solvers{1};
            %tic;
            [x1, fx1, exitflag1, output1] = package1(prob1);
            %T = toc; fprintf('\nRunning time for %s:\t %f\n', solvers{1}, T);
            prob2.options.solver = solvers{2};
            %tic;
            [x2, fx2, exitflag2, output2] = package2(prob2);
            %T = toc; fprintf('\nRunning time for %s:\t %f\n', solvers{2}, T);
        else
            prob1.options.solver = solvers{1};
            [x1, fx1, exitflag1, output1] = package1(prob1.objective, prob1.x0, prob1.Aineq, ...
                prob1.bineq, prob1.Aeq, prob1.beq, prob1.lb, prob1.ub, prob1.nonlcon, prob1.options);
            prob2.options.solver = solvers{2};
            [x2, fx2, exitflag2, output2] = package2(prob2.objective, prob2.x0, prob2.Aineq, ...
                prob2.bineq, prob2.Aeq, prob2.beq, prob2.lb, prob2.ub, prob2.nonlcon, prob2.options);
        end
    else
        if call_by_structure
            [x1, fx1, exitflag1, output1] = solver1(prob1);
            [x2, fx2, exitflag2, output2] = solver2(prob2);
        else
            switch lower(tested_solver_name)
            case {'uobyqa', 'newuoa'}
                [x1, fx1, exitflag1, output1] = solver1(prob1.objective, prob1.x0, prob1.options);
                [x2, fx2, exitflag2, output2] = solver2(prob2.objective, prob2.x0, prob2.options);
            case {'bobyqa'}
                [x1, fx1, exitflag1, output1] = solver1(prob1.objective, prob1.x0, prob1.lb, prob1.ub, prob1.options);
                [x2, fx2, exitflag2, output2] = solver2(prob2.objective, prob2.x0, prob2.lb, prob2.ub, prob2.options);
            case {'lincoa'}
                [x1, fx1, exitflag1, output1] = solver1(prob1.objective, prob1.x0, ...
                    prob1.Aineq, prob1.bineq, prob1.Aeq, prob1.beq, prob1.lb, prob1.ub, prob1.options);
                [x2, fx2, exitflag2, output2] = solver2(prob2.objective, prob2.x0, ...
                    prob2.Aineq, prob2.bineq, prob2.Aeq, prob2.beq, prob2.lb, prob2.ub, prob2.options);
            case {'cobyla'}
                [x1, fx1, exitflag1, output1] = solver1(prob1.objective, prob1.x0, ...
                    prob1.Aineq, prob1.bineq, prob1.Aeq, prob1.beq, prob1.lb, prob1.ub, prob1.nonlcon, prob1.options);
                [x2, fx2, exitflag2, output2] = solver2(prob2.objective, prob2.x0, ...
                    prob2.Aineq, prob2.bineq, prob2.Aeq, prob2.beq, prob2.lb, prob2.ub, prob2.nonlcon, prob2.options);
            otherwise
                error('Wrong solver tested: %s', tested_solver_name);
            end
        end
    end
catch exception
    % Do nothing for the moment
end

% Restore the random number generator state
rng(orig_rng_state);

if ~isempty(exception)
    if endsWith(exception.identifier, 'ConstraintFailureAtX0') && (strcmpi(solvers{1}, 'cobyla') || strcmpi(solvers{2}, 'cobyla'))
        % In this case, error is expected.
        equiv = true;
        return
    else
        rethrow(exception)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END: Call the solvers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

equiv = iseq(x1(:), fx1, exitflag1, output1, x2(:), fx2, exitflag2, output2, prec);

if ~equiv
    format long;
    fprintf('\nnf: nf1 = %d, nf2 = %d', output1.funcCount, output2.funcCount)
    fprintf('\nx:')
    x1(:)'
    x2(:)'
    (x1(:) == x2(:))'
    fprintf('\nf: fx1 = %.16e, fx2 = %.16e', fx1, fx2)
    fprintf('\nexitflag: exitflag1 = %d, exitflag2 = %d', exitflag1, exitflag2)
    nhist = min(length(output1.fhist), length(output2.fhist));
    fprintf('\nfhist (compare only the last %d evaluations):', nhist);
    output1.fhist
    output2.fhist
    fhist2 = output1.fhist(end-nhist+1: end);
    fhist1 = output2.fhist(end-nhist+1: end);
    fhist1 == fhist2
    if (isfield(output1, 'constrviolation'))
        fprintf('\nconstrviolation: constrviolation1 = %.16e, constrviolation2 = %.16e', ...
            output1.constrviolation, output2.constrviolation)
        fprintf('\nchist (compare only the last %d evaluations):', nhist);
        output1.chist
        output2.chist
        chist1 = output1.chist(end-nhist+1:end);
        chist2 = output2.chist(end-nhist+1:end);
        chist1 == chist2
    end
    if single_test && options.sequential
    %if options.sequential
        fprintf('\nThe solvers produce different results on %s at the %dth run.\n\n', pname, ir);
        cd(options.olddir);
        keyboard
    end
    error('\nThe solvers produce different results on %s at the %dth run.\n', pname, ir);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eq = iseq(x, f, exitflag, output, xx, ff, ee, oo, prec)
eq = true;

if ~isempty(setdiff(fieldnames(output), [fieldnames(oo); 'fhist'; 'xhist'; 'chist'; 'nlcihist'; 'nlcehist'])) ...
        || ~isempty(setdiff(fieldnames(oo), [fieldnames(output); 'fhist'; 'xhist'; 'chist'; 'nlcihist', 'nlcehist']))
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

if (norm(xx-x)/(1+norm(x)) > prec || abs(ff-f)/(1+abs(f)) > prec ...
        || abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation)) > prec)
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

minfhist = min(length(output.fhist), length(oo.fhist));
if norm(output.fhist(end-minfhist+1:end) - oo.fhist(end-minfhist+1:end))/(1+norm(output.fhist(end-minfhist+1:end))) > prec
    eq = false;
end

minchist = min(length(output.chist), length(oo.chist));
if norm(output.chist(end-minchist+1:end) - oo.chist(end-minchist+1:end))/(1+norm(output.chist(end-minchist+1:end))) > prec
    eq = false;
end

if (prec == 0 && (exitflag ~= ee|| oo.funcCount ~= output.funcCount))
    eq = false;
end

%diff = max([abs(ff-f)/(1+abs(f)), norm(xx-x)/(1+norm(x)), ...
%    abs(oo.constrviolation-output.constrviolation)/(1+abs(output.constrviolation))]);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function blacklist = black_list(solver)
%BLACK_LIST returns a list of problems that will be skipped when testing solver.
% Unless otherwise specified, the problems listed below take too much time and make the test
% on GitHub Actions run overtime.

blacklist = {};

% As of 20230426, the objective function of HS67 takes infinite time to be evaluated at some
% points, e.g., [88.1351318; 12829.9219; 1.0e-5], maybe due to an infinite cycling.
blacklist = [blacklist, {'HS67'}];

switch(lower(solver))
case 'uobyqa'
    blacklist = [blacklist, { ...
        'BENNETT5LS', ...
        'HATFLDFL', ...
        'HIELOW', ...
        'VARDIM', ...
        'YFITU', ...
        }];
case 'newuoa'
    blacklist = [blacklist, { ...
        'ARGLINB', ...
        'ARGLINC', ...
        'ARGTRIGLS', ...
        'BA-L1LS', ...
        'BA-L1SPLS', ...
        'BROWNAL', ...
        'COATING', ...
        'DIAMON3DLS', ...
        'DMN37143LS', ...
        'HYDC20LS', ...
        'HYDC20LS', ...
        'LUKSAN11LS', ...
        'LUKSAN13LS', ...
        'LUKSAN14LS', ...
        'LUKSAN15LS', ...
        'LUKSAN16LS', ...
        'LUKSAN17LS', ...
        'LUKSAN21LS', ...
        'LUKSAN22LS', ...
        'MANCINO', ...
        'PENALTY2', ...
        'PENALTY3', ...
        'QING', ...
        'VARDIM', ...
         }];
case 'bobyqa'
    blacklist = [blacklist, { ...
        'ARGLINA', ...
        'ARGLINB', ...
        'ARGLINC', ...
        'ARGTRIGLS', ...
        'BA-L1SPLS', ...
        'BROWNAL', ...
        'CHEBYQAD', ...
        'COATING', ...
        'DIAMON3DLS', ...
        'DMN15102LS', ...
        'DMN37143LS', ...
        'HOLMES', ...
        'LRA9A', ...
        'LUKSAN11LS', ...
        'LUKSAN12LS', ...
        'LUKSAN13LS', ...
        'LUKSAN14LS', ...
        'LUKSAN15LS', ...
        'LUKSAN16LS', ...
        'LUKSAN17LS', ...
        'LUKSAN21LS', ...
        'LUKSAN22LS', ...
        'LUKSAN23LS', ...
        'MANCINO', ...
        'PENALTY2', ...
        'PENALTY3', ...
        'QING', ...
        'VARDIM', ...
        }];
case 'lincoa'
    %blacklist = [blacklist, {'LSNNODOC', 'HS55', 'HEART6', 'AVGASA', 'AVGASB', 'CHEBYQAD', 'HS54'}]; % Classical lincoa encounters SEGFAULT
    blacklist = [blacklist, { ...
        'ARGTRIGLS', ...
        'BA-L1SPLS', ...
        'BQP1VAR', ...
        'BROWNAL', ...
        'CHEBYQAD', ...
        'COATING', ...
        'CVXQP1', ...
        'DECONVU', ...
        'DIAMON3DLS', ...
        'DMN15103LS', ...
        'DMN15332LS', ...
        'DMN37143LS', ...
        'DUAL1', ...
        'DUAL2', ...
        'DUAL3', ...
        'HIMMELBI', ...
        'HYDC20LS', ...
        'LINSPANH', ...
        'LRCOVTYPE', ...
        'LSQFIT', ...
        'LUKSAN11LS', ...
        'LUKSAN12LS', ...
        'LUKSAN13LS', ...
        'LUKSAN14LS', ...
        'LUKSAN15LS', ...
        'LUKSAN16LS', ...
        'LUKSAN17LS', ...
        'LUKSAN21LS', ...
        'LUKSAN22LS', ...
        'MANCINO', ...
        'MINSURF', ...
        'PENALTY3', ...
        'QING', ...
        'QPCBLEND', ...
        'QPCBOEI2', ...
        'QPNBLEND', ...
        'QPNBOEI2', ...
        'SIM2BQP', ...
        'SPANHYD', ...
        'VARDIM', ...
        }];
case 'cobyla'
    %blacklist = [blacklist, {'LAUNCH', 'MINMAXRB', 'MAKELA1', 'HS75', 'GAUSS3','HATFLDG'}]; % Classical cobyla encounters SEGFAULT
    %blacklist = [blacklist, {'EXTRASIM', 'POLAK2', 'POLAK6', 'SPIRAL'}]; % Assertion failed: B = A^{-1}
    blacklist = [blacklist, {'HS80'}];  % QRADD_RDIAG: Assertion failed: C^T*Q(:, N) == Rdiag(N).
    blacklist = [blacklist, {'DEGENLPA'}]; % Classical cobyla encounters infinite cycling
    blacklist = [blacklist, { ...
        'ACOPP30', ...
        'ACOPR14', ...
        'ACOPR30', ...
        'AIRPORT', ...
        'ANTWERP', ...
        'AVION2', ...
        'BATCH', ...
        'BQPGASIM', ...
        'BQPGABIM', ...
        'CHANDHEQ', ...
        'CHEBYQAD', ...
        'CHEBYQADNE', ...
        'CHNROSNB', ...
        'CHNRSBNE', ...
        'CHNRSNBMNE', ...
        'CORE1', ...
        'CRESC132', ...
        'DALLASS', ...
        'DECONVB', ...
        'DECONVBNE', ...
        'DECONVC', ...
        'DECONVU', ...
        'DEGENLPB', ...
        'DEGENQPC', ...
        'DIAMON2D', ...
        'DMN15102', ...
        'DMN15103', ...
        'DMN15332', ...
        'DMN15333', ...
        'DMN37142', ...
        'DMN37142LS', ...
        'DMN37143', ...
        'DNIEPER', ...
        'DUAL1', ...
        'DUAL2', ...
        'DUAL4', ...
        'DUALC5', ...
        'ERRINRSM', ...
        'ERRINRSMNE', ...
        'FBRAIN3', ...
        'FEEDLOC', ...
        'GROUPING', ...
        'HAIFAM', ...
        'HIMMELBI', ...
        'HS55', ...
        'HYDC20LS', ...
        'HYDCAR20', ...
        'KISSING2', ...
        'LAKES', ...
        'LINSPANH', ...
        'LOADBAL', ...
        'LUKSAN11', ...
        'LUKSAN11LS', ...
        'LUKSAN12', ...
        'LUKSAN12LS', ...
        'LUKSAN13', ...
        'LUKSAN13LS', ...
        'LUKSAN14', ...
        'LUKSAN14LS', ...
        'LUKSAN15', ...
        'LUKSAN16', ...
        'LUKSAN17', ...
        'LUKSAN17LS', ...
        'LUKSAN21', ...
        'LUKSAN21LS', ...
        'LUKSAN22', ...
        'LUKSAN22LS', ...
        'MANCINONE', ...
        'MESH', ...
        'METHANL8', ...
        'METHANL8LS', ...
        'MIFFLIN1', ...
        'MSS1', ...
        'NET1', ...
        'PALMER1NE', ...
        'PALMER4ANE', ...
        'PALMER5BNE', ...
        'PALMER7ANE', ...
        'PALMER8ENE', ...
        'POLAK2', ...
        'PRODPL0', ...
        'PRODPL1', ...
        'QPCBLEND', ...
        'SIPOW3', ...
        'SPANHYD', ...
        'SWOPF', ...
        'TAX13322', ...
        'TAXR13322', ...
        'TOINTQOR', ...
        'TRO4X4', ...
        'TRO6X2', ...
        'VANDERM1', ...
        'VANDERM2', ...
        'VANDERM3', ...
        'VANDERM4', ...
        'VESUVIOU', ...
        'WATER', ...
        }];
end
return
