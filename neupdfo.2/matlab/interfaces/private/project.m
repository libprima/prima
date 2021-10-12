function [x, fx, exitflag, output] = project (Ai, bi, Ae, be, lb, ub, x0, options)
%PROJECT is a (naive) semismooth-Newton method for solving 
%   min f(x) = 0.5*||x-x0||^2 
%   s.t. Ai*x <= bi
%   Ae*x = be
%   lb <= x <= ub
%   The method minimizes 
%   Fp(x) = f(x) + 0.5*sigma*(||Ae*x-be||^2+||(Ai*x-bi)_+||^2+||(lb-x)_+||^2+||(x-ub)_+||^2)  
%   by semismooth-Newton. Sigma is adaptively chosen. In theory, sigma
%   needs to tend to infinity in order to solve the problem (which is
%   a drawback of the Courant penalty function). 
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
%
% 1. Critial issues:
% 1.1. How to define sigma (intial value, maximal value, update)
% 1.2. The Semismooth-Newton method may have difficulty when the Jacobin
% is badly conditioned, which will happen when sigma is large.
%
% 2. exitflag:
% If output.algorithm = 'semismooth-newton':
% 0: optimal point found up to prescribed tolerance
% -1: non-optimal feasible point found 
% -2: no feasible point found (according to the prescribed tolerance, which may be very restrictive)
% If output.algorithm = 'quadprog' or 'fmincon', then the exitflag is the 
% exitflag of these MATLAB functions.
%
% 3. Although the code may be used to solve the above mentioned problem for
% general purposes, it was written to provide a feasilble starting point for
% Powell's LINCOA code, where the problem size is normally not big
% (typically tens of variables/constraints; at most thousands of them),
% and a feasible point is adequate. 
%
% 4. Testing results (without using 'quadprog' or 'fmincon'):
% For all the 235 bound/linearly constrained CUTEst problems with at 
% most 5,000 variables and 50,000 linear constraints, the code can find a 
% point with relative constration violation (RCV) at most 10^(-6) except the 
% following problems:
% a. BLOWEYA, BLOWEYB: final RCV between 10^(-5) and 10^(-6);
% b. DTOC3: final RCV between 10^(-5) and 10^(-4);
% c. POWELL20: final RCV between 10^(-4) and 10^(-3);
% d. HUES-MOD, HUESTIS: final RCV between 10^(-2) and 10^(-1);
% e. LINCONT: quadprog and fmincon of MATLAB cannot find a feasible point either
% f. ARGLALE, ARGLBLE, ARGLCLE, MODEL, NASH: problems infeasible accoring to quadprog
%
% 5. If no feasible point is found and TryMatlab=1, then we check whether the
% MATLAB on this machine has 'quadprog' or 'fmincon'. If yes, we use them to
% solve the problem. Be realistic. Our objective here is more to solve the
% problem than to develope a new algorithm! 
%
% TODO: Better algorithm/implemention for solving this projection problem 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% project starts

% Default options:
CTol = 1e-12; % Tolerance for relative constraint violation
OptTol = 1e-8; % Tolerance for first order optimality
FunTol = 1e-12; % Lower bound for the change in the value of Fp during a step
StepTol = 1e-10; % Lower bound for the length of a step
TryMatlab = 1;
Warnings = false;

maxit = 100;
mu = 0.4;
minalpha = 1e-15; % Lower bound for the step size of the Newton step
ainc = 1.5;
adec = 0.75;
cdec = 0.5;
cdecs = 0.05;
sinc = 1.2; % 1.2~1.5 seems good; 2 seems too big 
sdec = 0.6;
sigma = 10;
maxsigma = 1e15;
maxcon = 1e20; % This value will be used to decide whether an inequality constraint can be ignored or not.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin == 8 && isa(options, 'struct')) % Adopt user-defined options.
    if (isfield(options,'ConstraintTolerance'))
        CTol = options.ConstraintTolerance;
    end
    if (isfield(options,'OptimalityTolerance'))
        OptTol = options.OptimalityTolerance;
    end
    if (isfield(options,'FunctionTolerance'))
        FunTol = options.FunctionTolerance;
    end
    if (isfield(options,'StepTolerance'))
        StepTol = options.StepTolerance;
    end
    if (isfield(options,'MaxIterations')) 
        maxit = int32(options.MaxIterations);
    end
    if (isfield(options,'TryMatlab')) 
        TryMatlab = options.TryMatlab;
    end
    if (isfield(options,'Warnings')) 
        Warnings = options.Warnings;
    end
end

if ~Warnings % Suppress the warning about nearly singular matrices
    orig_warning_state = warning;
    warning('off','MATLAB:nearlySingularMatrix');
end

n = length(x0);
if(isempty(lb))
    lb = -inf(n,1);
end
if(isempty(ub))
    ub = inf(n,1);
end

if((isempty(Ai) || min(bi) >= norm(Ai, inf)*maxcon) && max(lb) <= -maxcon && min(ub) >= maxcon) % There are only equality constraints
    if(~isempty(Ae))
        if ~isempty(which('lsqminnorm'))
            x = lsqminnorm(Ae, be-Ae*x0) + x0;
        else
            x = Ae\(be-Ae*x0) + x0; 
        end
        x = min(max(x, lb), ub);
        fx = 0.5*(x-x0)'*(x-x0);
        output.penalty = 0;
        output.iterations = 0; 
        if(isempty(Ae))
            re = 0;
            be = 0;
        else
            re = Ae*x-be;
        end
        if(isempty(Ai))
            ri = 0;
            bi = 0;
        else
            ri = max(Ai*x-bi, 0);
        end
        output.constrviolation = norm([re; ri])/(1+ norm([be; bi])); 
        if (output.constrviolation > CTol)
            exitflag = -2;
        else
            exitflag = 0;
        end
        return;
    end
end

x = min(max(x0, lb), ub);
if((isempty(Ai) || min(bi) >= norm(Ai, inf)*maxcon) && isempty(Ae)) % There are only bounds
    fx = 0.5*(x-x0)'*(x-x0);
    output.penalty = 0;
    output.iterations = 0; 
    if(~isempty(Ai))
        output.constrviolation = norm([0; max(Ai*x-bi,0)])/(1+norm([0; bi]));
    else
        output.constrviolation = 0;
    end
    if (output.constrviolation > CTol)
        exitflag = -2;
    else
        exitflag = 0;
    end
    return;
end

if(isempty(Ae))
    Ae = zeros(1,n);
    be = 0;
end
if(isempty(Ai))
    Ai = zeros(1,n);
    bi = 0;
end

I = speye(n,n);
b = -x0;
cv = NaN;
befin = [0; be(abs(be)<inf)]; % Put a zero to avoid empty array
bifin = [0; bi(abs(bi)<inf)]; 
ubfin = [0; ub(abs(ub)<inf)];
lbfin = [0; lb(abs(lb)<inf)];
denomc = 1 + norm(befin) + norm(bifin) + norm(ubfin) + norm(lbfin);

for k = 1 : maxit
    re = Ae*x-be;
    ri = Ai*x-bi;
    rlb = -x+lb;
    rub = x-ub;
    cv_old = cv;
    cv = norm([re; max(ri, 0); max(rlb, 0); max(rub, 0)]);
    if (cv > max(0.1*CTol*denomc, cdec*cv_old))
        sigma = min(sinc*sigma, maxsigma);
    elseif (cv < cdecs*cv_old)
        sigma = sdec*sigma;
    end
    g = x+b;
    ge = Ae'*re;
    gi = Ai'*max(ri,0);
    glb = - max(rlb,0);
    gub = max(rub,0);
    Pg = g + sigma*(ge + gi + glb + gub);
    denomg = 1+norm(g) + sigma*(norm(ge+gi+glb+gub));
    if (norm(Pg)/denomg <= OptTol)
        if (cv/denomc <= CTol)
            break;
        else
            continue;
        end
    end
    Ji = Ai(ri>0, :)'* Ai(ri>0, :);
    if(isempty(Ji))
        J = I + sigma*(Ae'*Ae + double(diag(rlb>0)) + double(diag(rub>0)));
    else
        J = I + sigma*(Ae'*Ae + Ji + double(diag(rlb>0)) + double(diag(rub>0)));
    end
    % There should be much better ways to solve the linear equation J*d = -g !
    % We are lazy here. 
    d = -J\Pg;
    if (norm(J*d+Pg)/norm(Pg) >= 1e-1) && ~isempty(which('lsqminnorm'))
        d = lsqminnorm(-J, Pg);
    end
    dPg = d'*Pg;
    F = @(x) Fp(x, sigma, b, Ai, bi, Ae, be, lb, ub);
    Fx = F(x);
    if (dPg >= -FunTol*max(1, abs(Fx)))
        continue;
    end
    normx = max(norm(x), 1);
    normd = norm(d);
    alpha = 1;
    if (F(x+alpha*d)- Fx <= mu*alpha*dPg)
        while (F(x+ainc*alpha*d) - Fx <= mu*ainc*alpha*dPg) 
            alpha = ainc*alpha;
        end
    else
        while (F(x+alpha*d) - Fx > mu*alpha*dPg && alpha > minalpha && alpha*normd > StepTol*normx && -alpha*dPg > FunTol*max(abs(Fx),1)) 
            alpha = adec*alpha;
        end
    end
    if (alpha > minalpha && alpha*normd > StepTol*normx && -alpha*dPg > FunTol*max(abs(Fx),1))
        x = x + alpha*d;
    elseif (sigma == maxsigma)
        break;            
    end
end

re = Ae*x-be;
ri = Ai*x-bi;
rlb = -x+lb;
rub = x-ub;
cv = norm([re; max(ri, 0); max(rlb, 0); max(rub, 0)]);
Pg = g + sigma*(Ae'*re + Ai'*max(ri,0) - max(rlb,0) + max(rub,0));
denomg = 1 + norm(x+b) + sigma*(norm(Ae'*(Ae*x-be)+Ai'*(max(Ai*x-bi, 0))-max(lb-x, 0)+max(x-ub, 0)));
if (cv > CTol*denomc)
    exitflag = -2;
elseif (norm(Pg) > OptTol*denomg)
    exitflag = -1;
else
    exitflag = 0;
end
fx = 0.5*(x-x0)'*(x-x0);
output.algorithm = 'semismooth-newton';
output.penalty = sigma; 
output.iterations = k; 
output.firstorderopt = norm(Pg)/denomg + cv/denomc;
output.constrviolation = cv/denomc;

x_save = x;
fx_save = fx;
output_save = output;
exitflag_save = exitflag;

MatlabVersion = ver;
TryMatlab = TryMatlab && any(strcmp({MatlabVersion.Name}, 'Optimization Toolbox')); % Do not try matlab unless Optimization Toolbox is available
if output.constrviolation > 10*CTol && TryMatlab % No feasible ponit was found. Try quadprog or fmincon. 
    if ~isempty(which('quadprog'))
        options = optimoptions('quadprog'); 
        options.Display = 'off'; % No talking 
        options.ConstraintTolerance = CTol;
        options.MaxIterations = maxit;
%        [~, x, fx, exitflag, output] = evalc('quadprog (speye(n,n), -x0, Ai, bi, Ae, be, lb, ub, x0, options)'); % We do not want any message printed by quadprog
        [x, fx, exitflag, output] = quadprog (speye(n,n), -x0, Ai, bi, Ae, be, lb, ub, x0, options); % quadprog is silent with the Display set to 'off'
        fx = fx + 0.5*(x0'*x0);
        output.algorithm = 'quadprog (MATLAB)';
    elseif ~isempty(which('fmincon'))
        options = optimoptions('fmincon'); 
        options.Display = 'off'; % No talking
        options.SpecifyObjectiveGradient = true;
        options.Algorithm = 'sqp';
        options.ConstraintTolerance = CTol;
        options.MaxIterations = maxit;
%        [~, x, ~, exitflag, output] = evalc('fmincon(@(x) dist_sq(x, x0), x0, Ai, bi, Ae, be, lb, ub, [], options)'); % We do not want any message printed by fmincon   
%        fx = dist_sq(x, x0); % We may return fx when calling fmincon. However, without this line, checkcode will complain that dist_sq is unused due to evalc. 
        [x, fx, exitflag, output] = fmincon(@(x) dist_sq(x, x0), x0, Ai, bi, Ae, be, lb, ub, [], options); % fmincon is silent with the Display set to 'off'
        output.algorithm = 'fmincon (MATLAB)';
    end
end

if(isempty(x)) % x can be [] if quadprog thinks the problem is infeasible (which may not be true; for example, HIMMELBJ problem in CUTEst)
    x = x_save;
    fx = fx_save;
    output = output_save;
    exitflag = exitflag_save;
end

if ~Warnings
    warning(orig_warning_state); % Restore the behavior of displaying warnings
end

% project ends
return

function Fval = Fp(x, sigma, b, Ai, bi, Ae, be, lb, ub) 
re = Ae*x-be;
ri = Ai*x-bi;
rlb = -x+lb;
rub = x-ub;
rip = [0; ri(ri>0)]; % Put a zero to avoid empty array
rlbp = [0; rlb(rlb>0)];
rubp = [0; rub(rub>0)];
Fval = (0.5*(x'*x)+b'*x) + 0.5*sigma*(re'*re + rip'*rip + rlbp'*rlbp + rubp'*rubp);
return

function [f, g] = dist_sq(x, x0)
f = 0.5*(x-x0)'*(x-x0);
g = x-x0;
return
