function frec = testalg(solver, problem, n, options)

% This version is not intended to be released. It is only for test. 
%
% All rights reserved. 
%
% ZHANG Zaikun, 08/08/2016
% Department of Applied Mathematics, The Hong Kong Polytechnic University

global fval_history; 

if (nargin == 4 && isfield(options, 'reproduce') && options.reproduce == true) % Set options.reproduce = true for reproducing results of the last experiment. 
    load('seed.testalg.mat', 'seed');
else
    rng('shuffle'); 
    seed = 1e8*(2*rand - 1); 
    save('seed.testalg.mat', 'seed');
end

% Remark on rng: all the calls of rng are to make sure that
% 1. The randomness is the same for different solvers.
% 2. The randomness is independent of each other for differen problems (indexed by ip).
% 3. The randomness is independent of each other for different runs (indexed by ir).
% 4. The randomness is independent of each other for different calls of testalg.
% 5. The randomness is nearly independent of each other for different x.
% 6. The randomness is reproducible. 

rhoend = 1e-6;
tol = rhoend;
maxfun = 500*n;
maxiter = maxfun; 
ftarget = -Inf;
randomizex0 = 0;
noiselevel = 0;
dnoiselevel = 0;
nr = 5;

if (nargin == 4)
    if (isfield(options, 'tol'))
        tol = options.tol;
        rhoend = tol;
    end
    if (isfield(options, 'rhoend'))
        rhoend = options.rhoend;
        tol = rhoend;     
    end
    if (isfield(options, 'maxfun'))
        maxfun = options.maxfun;
    end
    if (isfield(options, 'maxiter'))
        maxiter = options.maxiter;
    end
    if (isfield(options, 'ftarget'))
        ftarget = options.ftarget;
    end
    if (isfield(options, 'randomizex0'))
        randomizex0 = options.randomizex0;
    end
    if (isfield(options, 'randrun'))
        nr = options.randrun;
    end
    if (isfield(options, 'noise'))
        if (isnumeric(options.noise)) 
            noi.type = 'relative';
            noi.level = options.noise;
            options = rmfield(options, 'noise');
            options.noise = noi;
        end
        noiselevel = options.noise.level;
    end
    if (isfield(options, 'dnoise'))
        if (isnumeric(options.dnoise)) 
            dnoi.type = 'relative';
            dnoi.level = options.dnoise;
            options = rmfield(options, 'dnoise');
            options.dnoise = dnoi;
        end
        dnoiselevel = options.dnoise.level;
    end
else
    options = [];
end

if (randomizex0 == 0 && noiselevel == 0)
    nr = 1;
end 

if (strcmpi('ALL', problem))
    prob = textread('problems', '%s');
else
    prob = {problem};
end
np = length(prob);

if (strcmpi('ALL', solver))
    sol = textread('solvers', '%s');
else
    sol = {solver};
end
ns = length(sol);

frec = NaN(np, ns, nr, maxfun);
for ip = 1 : np 
    display('+++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    display(strcat(int2str(ip), '. ', prob{ip}, ':'));
    [x0, rhobeg] = setuptest(prob{ip}, n);
    for is = 1 : ns 
        solv = sol{is};
        if (np*nr <= 20)
            display('-------------------------------------------------------');
        end
        display(strcat(sol{is}, ':'));
        for ir = 1 : nr 
            rng(int32(1e8*abs(sin(seed) + sin(double(ir)) + sin(sum(double(prob{ip}))) + sin(double(n)) + sin(1e8*randomizex0))));
            if (strcmpi('chebquad', prob{ip}))
                xr = x0 + min(abs(randomizex0), 1e-1)*(2*rand(n,1) - 1).*max(ones(n, 1), abs(x0)); 
            else
                xr = x0 + randomizex0*randn(n,1).*max(ones(n, 1), abs(x0)); 
            end
            if (exist('solv_options', 'var') == 1)  
                clear solv_options;
            end
            solv_options.rhobeg = rhobeg;
            solv_options.rhoend = rhoend;
            solv_options.tol = tol;
            solv_options.maxfun = maxfun;
            solv_options.maxiter = maxiter;
            solv_options.maxit = maxiter; 
            solv_options.ftarget = ftarget;
            solv_options.npt = 2*n+1;
            if (isfield(options, 'debug')) 
                solv_options.debug = options.debug;
            end
            if (isfield(options, 'chkfunval'))
                solv_options.chkfunval = options.chkfunval;
            end
            solv_options.warndim = false;
            if (isfield(options, 'reproduce'))
                solv_options.reproduce = options.reproduce;
            end
            if (isfield(options, 'signif')) 
                solv_options.signif = options.signif; 
            end
            if (isfield(options, 'prec')) 
                solv_options.prec = options.prec; 
            end
            if (isfield(options, 'noise')) 
                solv_options.noise = options.noise; 
            end
            if (isfield(options, 'dnoise')) 
                solv_options.dnoise = options.dnoise; 
            end

            if (strcmpi(solv, 'newuoasdev3'))
                solv = 'newuoasdev';
                solv_options.maxsubspacedim = 3;
            elseif (strcmpi(solv, 'newuoasdev4'))
                solv = 'newuoasdev';
                solv_options.maxsubspacedim = min(n,4);
            elseif (strcmpi(solv, 'newuoasdev5'))
                solv = 'newuoasdev';
                solv_options.maxsubspacedim = min(n,5);
            elseif (strcmpi(solv, 'newuoasdev6'))
                solv = 'newuoasdev';
                solv_options.maxsubspacedim = min(n,6);
            elseif (strcmpi(solv, 'newuoasdev10'))
                solv = 'newuoasdev';
                solv_options.maxsubspacedim = min(n,10);
            elseif (strcmpi(solv, 'newuoasdev20'))
                solv = 'newuoasdev';
                solv_options.maxsubspacedim = min(n,20);
            end

            if (exist('fun_options', 'var') == 1)  
                clear fun_options;
            end
            fun_options = options;
            if (isfield(fun_options, 'noise'))
                rng(int32(1e8*abs(sin(seed) + sin(double(ir)))));
                fun_options.noise.seed = 1e8*(2*rand - 1);
            end

            fval_history = []; 
            x = feval(solv, @(x)evalobjfun(prob{ip}, x, fun_options), xr, solv_options);
            nf = length(fval_history);
            nf = min(nf, maxfun); % Some solvers do not respect maxfun strictly, for example, fminunc of matlab.
            frec(ip, is, ir, 1:nf) = fval_history(1:nf);
            frec(ip, is, ir, nf+1:maxfun) = fval_history(nf);
            if (np*nr <= 20)
                fprintf('nf = %d\n', nf);
                minfval = min(fval_history);
                fprintf('minfval = %.6E\n', minfval);
                if (n <= 30)
                    xx=sprintf('%.3E ', x');
                    fprintf('xbest = %s\n', xx);
                end
            end
        end
    end
end

return;


function f = evalobjfun(fun, x, options)

%options:
%regul: regul.lambda regul.p 
%noise: noise.type noise.level noise.seed
%prec
%signif

global fval_history;

% Calculate the accurate function value, and save it to fval_history.
f = evalfun(fun, x); 
if (isfield(options, 'regul'))
    if (isfield(options.regul, 'lambda'))
        lambda = options.regul.lambda;
        p = options.regul.p;
    else
        lambda = 10;
        p = options.regul;
    end
    if (p >= 0) 
        f = f + lambda*sum(abs(x).^p);
    end
end
fval_history = [fval_history, f];

% Now calculate the function value that will be fed to the solvers.
if (isfield(options, 'noise'))
    noise = options.noise;
    if (abs(noise.level) > 0)
        rng(int32(1e8*abs(sin(double(noise.seed)) + sin(1e8*noise.level) + sin(sum(abs(sin(1e8*x)))) + sin(double(length(x))) + sin(sum(double(fun))))));
        if (~isfield(noise, 'type'))
            noise.type = 'relative';
        end
        if (strcmpi(noise.type, 'absolute') || strcmpi(noise.type, 'additive') || strcmpi(noise.type, 'add') || strcmpi(noise.type, 'a') || strcmpi(noise.type, '+'))
            f = f + noise.level*randn;
        else
            f = f * (1 + noise.level*randn);
        end
    end
end

if (isfield(options, 'dnoise'))
    dnoise = options.dnoise;
    if (abs(dnoise.level) > 0)
        %phi0 = 0.6*cos(1e8*norm(x,9)) + 0.3*sin(100*norm(x,1))*cos(100*norm(x,Inf)) + 0.1*cos(norm(x));
        phi0 = 0.9*sin(100*norm(x,1))*cos(100*norm(x,Inf)) + 0.1*cos(norm(x));
        noisimul = phi0*(4*phi0^2-3);
        if (~isfield(dnoise, 'type'))
            dnoise.type = 'relative';
        end
        if (strcmpi(dnoise.type, 'absolute') || strcmpi(dnoise.type, 'additive') || strcmpi(dnoise.type, 'add') || strcmpi(dnoise.type, 'a') || strcmpi(dnoise.type, '+'))
            f = f + dnoise.level*noisimul;
        else
            f = f * (1 + dnoise.level*noisimul);
        end
    end
end

if (isfield(options, 'prec'))
    if (strcmpi(options.prec, 'single') || strcmpi(options.prec, 's'))
        x = double(single(x));
        f = single(evalfun(fun, x)); 
        if (isfield(options, 'regul'))
            if (isfield(options.regul, 'lambda'))
                lambda = options.regul.lambda;
                p = options.regul.p;
            else
                lambda = 10;
                p = options.regul;
            end
            if (p >= 0) 
                f = f + single(lambda*sum(abs(x).^p));
            end
        end
        f = double(f);
    end
end

if (isfield(options, 'signif'))
    sig = min(max(1, int32(options.signif)), 16); 
    %sf = round(f, sig, 'significant');
    sf = eval(mat2str(f, sig)); 
    f = sf + (f-sf)*(1-sin(sin(double(sig)) + sin(1e8*f) + sum(abs(sin(1e8*x)) + sin(double(length(x))) + sin(sum(double(fun)))))); % This makes the truncation more "irregular". 
end

return;
