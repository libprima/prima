function testpdfo(release, precision, nrun) 
%TESTPDFO tests pdfo on a few VERY simple problems.
%
%   Note: Do NOT follow the syntax here when you use pdfo. This file is 
%   written for testing purpose, and it uses quite atypical syntax. See 
%   rosenbrock_example.m for an illustration about how to use pdfo. 
%
%   ***********************************************************************
%   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk) 
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%
%   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
%   ***********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attribute: public (can be called directly by users)
% 
% TODO: None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_solvers =  {'uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla', 'pdfo'};

% Supress the following warning
orig_warning_state = warning;
cellfun(@(solver) warning('off', [solver, ':Debug']), all_solvers);
cellfun(@(solver) warning('off', [solver, ':ChkFunval']), all_solvers);
cellfun(@(solver) warning('off', [solver, ':Classical']), all_solvers);
cellfun(@(solver) warning('off', [solver, ':ReviseX0']), all_solvers);

if nargin < 1
    release = true;
end
if nargin < 2
    precision = 1e-6; % Default testing precision
end
if nargin < 3
    nrun = 1; % Number of runs with randomly perturbed x0
    perturb = 0; % Magnitude of permission on x0
else
    perturb = 1e-15;
end

options.debug = true;
options.chkfunval = true;

type_list = {'unconstrained', 'bound-constrained', 'linearly-constrained', 'nonlinearly-constrained'};

solver_list = {{'uobyqa', 'newuoa', 'bobyqa', 'lincoa', 'cobyla'}, ...
               {'bobyqa', 'lincoa', 'cobyla'}, ...
               {'lincoa', 'cobyla'}, ...
               {'cobyla'}
              };

fun_list = {@goldp, @mcc, @hmlb, @chrosen, @chebquad};

x0_list = {zeros(2,1), zeros(2,1), zeros(2,1), zeros(3,1), zeros(4,1)};

fopt_list = {{3, 3.2577598292927021e1, 600, 600,}, ... % goldp
             {-1.913222954981037, -5.0768934966810497e-01, 0.934627336363561,  0.934627336363561}, ... %mcc
             {0, 144.125, 136, 1.293402195198189e+02}, ... %hmlb
             {0, 0.7, 6.0153819489865101e-01, 2.4054273877510723e-01}, ... %chrosen
			 {0, 2.1272916538809911e-01,  2.9733407628744646e-01, 2.6143791089220641e-02} %chebquad
			};

%xopt_list = {{[0;-1], [-0.5; -4.9018111037601469e-01], [0; 0], [0; 0]}, ... % goldp 
%             {[-0.547197554599523; -1.547197552199337], [-1.4849642130761465e-01; -0.5], [0.267825716190182; 0], [0.267825716190182; 0]}, ... %mcc
%			  {[3; 2], [0.5; 0.5], [1; 0], [0.773622934941950; 0.633646237684564]}, ... %hmlb
%             {ones(3,1), ...
%             [0.5; 0.4; 0.16], ...
%             [6.0367489816330211e-01; 3.6090702088918725e-01; 3.5418080947510699e-02], ...
%             [7.5499486749647116e-01; 5.8016552326623394e-01; 3.0559894579354657e-01]}, ... %chrosen
%             {[5.9379627877490138e-01; 4.0620384108718871e-01; 8.9732723150433336e-01; 1.0267283287701645e-01], ...
%			  [1.0020951281679348e-01; 0.5; 1.0020960717311733e-01; 0.5], ... 
%			  [6.3140858383707379e-02; 6.3140842748373699e-02; 3.0073419977477961e-01; 5.7298409909313941e-01], ... 
%			  [3.2005844354004209e-01; 6.6425556098486729e-02; 4.9430560542604501e-01; 8.0548879983696642e-01]} %chebquad
%            }; 
% Note: 
% 1. chebquad is invariant with respect to permutations of the variables. Thus there is no unique xopt.
% 2. Himmelblau's function (hmlb) has multiple minima when
% unconstrained. They are [3; 2], [-2.805118; 3.131312], [-3.779310; -3.283186], [3.584428; -1.848126]

if release
    clflag_list = {false}; % clflag: call the solvers in classical mode (true) or not (false)
else
    clflag_list = {true, false}; % clflag: call the solvers in classical mode (true) or not (false)
end

pass = true; % All tests are successful or not

for irun = 1 : nrun
    fprintf ('\n');
    if (nrun > 1)
        fprintf ('Test %d:\n\n', irun);
    end
    try
        for itype = 1 : length(type_list)
            type = type_list{itype};
            fprintf ('Testing %s problems ...\n', strrep(type, '-', ' '));
            success = true; % The current test is successful or not
            for isolver = 1 : length(solver_list{itype})
                solver = solver_list{itype}{isolver};
                for ifun = 1 : length(fun_list)
                    for iclflag = 1 : length(clflag_list)
                        clflag = clflag_list{iclflag};
                        fun = fun_list{ifun};
                        x0 = x0_list{ifun};
                        if ~strcmp(solver, 'cobyla') || ~strcmp(func2str(fun), 'chebquad') 
                            % The result of cobyla on chebquad is sensitive to x0, so we do not perturb x0
                            r = abs(sin(1e3*sum(double([solver, func2str(fun), type]))*irun*(1:length(x0))'));
                            % Introduce a tiny perturbation to the experiments.
                            % We use a deterministic permutation so that
                            % experiments can be easily repeated when necessary.
                            x0 = x0 + perturb*r;
                        end
                        %xopt = xopt_list{ifun}{itype};
    				    fopt = fopt_list{ifun}{itype};
                        n = length(x0);
    
                        problem = struct();
                        problem.objective = fun;
                        problem.x0 = x0;
                        options.solver = solver;
                        options.classical = clflag;
                        problem.options = options;
    
                        switch type
                        case 'unconstrained'
                            [x, fx] = pdfo(fun, x0, options);
                        case 'bound-constrained' % 0 <= x <= 0.5
                            lb = -0.5*ones(n,1);
                            ub = 0.5*ones(n,1);
                            [x, fx] = pdfo(fun, x0, [], [], [], [], lb, ub, options);
                            problem.lb = lb;
                            problem.ub = ub;
                        case 'linearly-constrained' % simplex 
                            Aineq = ones(1, n);
                            bineq = 1;
                            lb = zeros(n,1);
                            [x, fx] = pdfo(fun, x0, Aineq, bineq, [], [], lb, [], options); 
                            problem.Aineq = Aineq;
                            problem.bineq = bineq;
                            problem.lb = lb;
                        case 'nonlinearly-constrained' % intersection of the unit ball and the positive orthant
                            nonlcon = @(x) ballcon(x, zeros(n,1), 1);
                            lb = zeros(n,1);
                            [x, fx] = pdfo(fun, x0, [], [], [], [], lb, [], nonlcon, options);
                            problem.lb = lb;
                            problem.nonlcon = nonlcon;
                        end
    
                        xs = pdfo(problem);
    
                        if strcmp(solver, 'cobyla') % The precision of cobyla is lower
                            prec = max(1e3*precision, 1e-2);
                        else
                            prec = precision;
                        end
                        if ~release 
                            fprintf('\nsolver = %s,\tfun = %s,\t\tfx = %.16e,\t\tfopt = %.16e', solver, func2str(fun), fx, fopt);
                        end
                        if (norm(x-xs) > 0) || ((fx-fopt)/max(1, abs(fopt)) > prec) || (~release && abs(fx-fopt)/max(1, abs(fopt)) > prec)
                            success = false;
                            if clflag  
                                fprintf('\n**** pdfo (classical mode) FAILED a test on %s problem: solver = ''%s'', objective function = ''%s''.\n', type, solver, func2str(fun));
                            else
                                fprintf('\n**** pdfo FAILED a test on %s problem: solver = ''%s'', objective function = ''%s''.\n', type, solver, func2str(fun));
                            end
                            fprintf ('Required precision = %.2e,\t\tactual precision = %.2e\n', prec, abs(fx-fopt)/max(1, abs(fopt)));
                        end
                    end
                end
            end
            if success 
                if ~release
                    fprintf('\n\n');
                end
                fprintf ('Succeed.\n\n');
    		else
            	pass = false;
                fprintf ('\n');
            end
        end
    catch exception
        fprintf('\nERROR occurred during the test:\n\n');
        rethrow(exception);
    end
    
    if pass 
        fprintf('All tests were successful.\n\n');
    else    
        fprintf('Some tests FAILED.\n\n');
        break;
    end
end

warning(orig_warning_state); % Restore the behavior of displaying warnings

return

function [f, g, H]=chrosen(x)
%CHROSEN calculates the function value, gradient, and Hessian of the
%   Chained Rosenbrock function.
%   See 
%   [1] Toint (1978), 'Some numerical results using a sparse matrix
%   updating formula in unconstrained optimization'
%   [2] Powell (2006), 'The NEWUOA software for unconstrained
%   optimization without derivatives'

n=length(x);

alpha = 4;

f=0; % Function value
g=zeros(n,1); % Gradient
H=zeros(n,n); % Hessian

for i=1:n-1
    f = f + (x(i)-1)^2+alpha*(x(i)^2-x(i+1))^2;
  
    g(i)   = g(i) + 2*(x(i)-1)+alpha*2*(x(i)^2-x(i+1))*2*x(i);
    g(i+1) = g(i+1) - alpha*2*(x(i)^2-x(i+1));
  
    H(i,i)    =  H(i,i)+2+alpha*2*2*(3*x(i)^2-x(i+1));
    H(i,i+1)  =  H(i,i+1)-alpha*2*2*x(i);
    H(i+1,i)  =  H(i+1,i) -alpha*2*2*x(i);
    H(i+1,i+1)=  H(i+1,i+1)+alpha*2;
end

return

function f = chebquad(x)
%CHEBQUAD evaluates the Chebyquad function.
%
%   See 
%   [1] Fletcher (1965), 'Function minimization without evaluating derivatives --- a review'

n = length(x);
y(1,1:n) = 1;
y(2, 1:n) = 2*x(1:n) - 1;
for i = 2:n
    y(i+1, 1:n) = 2*y(2, 1:n).*y(i, 1:n) - y(i-1, 1:n);
end
f = 0;
for i = 1 : n+1
    tmp = mean(y(i, 1:n));
    if (mod(i, 2) == 1) 
        tmp=tmp+1/double(i*i-2*i);
    end
    f = f + tmp*tmp;
end

return

function [f, g] = hmlb(x)
%HMLB evaluates the Himmelblau's function and its gradient
%
%   See 
%   [1]  Himmelblau (1972),  'Applied Nonlinear Programming'

f = (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;
g = 2*[-7 + x(1) + x(2)^2 + 2*x(1)*(-11 + x(1)^2 + x(2)); -11 + x(1)^2 + x(2) + 2*x(2)*(-7 + x(1) + x(2)^2)];

return

function f = goldp(x)
%GOLDP evaluates the Goldstein-Price function
%
%   See
%   [1] Dixon, L. C. W., & Szego, G. P. (1978). The global optimization problem: an introduction. Towards global optimization, 2, 1-15.

f1a = (x(1) + x(2) + 1)^2;
f1b = 19 - 14*x(1) + 3*x(1)^2 - 14*x(2) + 6*x(1)*x(2) + 3*x(2)^2;
f1 = 1 + f1a*f1b;

f2a = (2*x(1) - 3*x(2))^2;
f2b = 18 - 32*x(1) + 12*x(1)^2 + 48*x(2) - 36*x(1)*x(2) + 27*x(2)^2;
f2 = 30 + f2a*f2b;

f = f1*f2;

return

function f = mcc(x)
%MCC evaluates the McCormick function

f1 = sin(x(1) + x(2));
f2 = (x(1) - x(2))^2;
f3 = -1.5*x(1);
f4 = 2.5*x(2);

f = f1 + f2 + f3 + f4 + 1;

return

function [cineq, ceq] = ballcon(x, centre, radius)
% BALLCON represnts the ball constraint ||x-centre|| <= radius

cineq = (x-centre)'*(x-centre) - radius^2; 
ceq = [];

return
