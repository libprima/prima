function tough_problem = tough(problem, random_seed, noise_level, with_failure)
%This function prepares the TOUGH version of a given problem.
% problem: a structure defining the original problem
% random_seed: a seed provided by the caller in order to ensure reproducibility
% noise_level: level of the noise
% with_failure: whether to fail the objective and constraint evaluation randomly

if nargin < 3
    noise_level = 2.0e-1;  % The noise level.
end
if nargin < 4
    with_failure = true;  % Whether to fail the function evaluation randomly.
end

% Set the random seed
orig_rng_state = rng();
rng(random_seed);

% Copy the problem options
if isfield(problem, 'options')
    tough_problem.options = problem.options;
else
    tough_problem.options = [];
end

% Set the starting point
x0 = problem.x0;
n = length(x0);
tough_problem.x0 = x0 + noise_level * max(1, abs(x0)) .* randn(n,1);

% Set the objective function
tough_problem.objective = @(x) tough_feval(problem.objective, x, random_seed, noise_level, with_failure);

% Set the bound constraints
minlb = -1.0e10;
maxub = 1.0e10;
lb = problem.lb;
ub = problem.ub;
if isempty(lb)
    lb = zeros(n, 1) + minlb;
else
    lb = max(minlb, lb);
end
if isempty(ub)
    ub = zeros(n, 1) + maxub;
else
    ub = min(maxub, ub);
end
if isempty(problem.lb)
    tough_problem.lb = [];
else
    tough_problem.lb = lb + noise_level * 0.9 * min(1, (ub - lb)) .* (rand(n, 1) - 1);
end
if isempty(problem.ub)
    tough_problem.ub = [];
else
    tough_problem.ub = ub + noise_level * 0.9 * min(1, (ub - lb)) .* (rand(n, 1) - 1);
end

% Set the linear constraints
Aeq = problem.Aeq;
beq = problem.beq;
if isempty(Aeq)
    tough_problem.Aeq = zeros(0, n);
    tough_problem.beq = zeros(0, 1);
else
    tough_problem.Aeq = Aeq + noise_level * max(1, abs(Aeq)) .* randn(size(Aeq, 1), size(Aeq, 2));
    tough_problem.beq = beq + noise_level * max(1, abs(beq)) .* randn(length(beq), 1);
end
Aineq = problem.Aineq;
bineq = problem.bineq;
if isempty(Aineq)
    tough_problem.Aineq = zeros(0, n);
    tough_problem.bineq = zeros(0, 1);
else
    tough_problem.Aineq = Aineq + noise_level * max(1, abs(Aineq)) .* randn(size(Aineq, 1), size(Aineq, 2));
    tough_problem.bineq = bineq + noise_level * max(1, abs(bineq)) .* randn(length(bineq), 1);
end

% Set the nonlinear constraints
if isempty(problem.nonlcon)
    tough_problem.nonlcon = [];
else
    tough_problem.nonlcon = @(x) tough_ceval(problem.nonlcon, x, random_seed, noise_level, with_failure);
end

% Restore the random seed
rng(orig_rng_state);

% `tough` ends here
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = tough_feval(func, x, random_seed, noise_level, with_failure)
%This function evaluates the function func at x for the TOUGH test.
if nargin < 4
    noise_level = 2e-1;
end
if nargin < 5
    with_failure = true;
end
f = func(x);
f = contaminate(f, x, random_seed, noise_level, with_failure);

% `tough_feval` ends here
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cineq, ceq] = tough_ceval(con, x, random_seed, noise_level, with_failure)
%This function evaluates the constraint function con at x for the TOUGH test. Here, following the
% convention of MATLAB, the constraint function is assumed to be of the form [cineq, ceq] = con(x),
% and the constraint is cineq <= 0 and ceq == 0.
if nargin < 4
    noise_level = 2e-1;
end
if nargin < 5
    with_failure = true;
end
[cineq, ceq] = con(x);
cineq = arrayfun(@(c) contaminate(c, x, random_seed, noise_level, with_failure), cineq);
ceq = arrayfun(@(c) contaminate(c, x, random_seed, noise_level, with_failure), ceq);

% `tough_ceval` ends here
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = contaminate(f, x, random_seed, noise_level, with_failure)
%This function contaminates f for the TOUGH test.
% f is the function value to be contaminated.
% x is the value of the decision variable corresponding to f.
% The random seed used internally (see `rseed` below) will be defined by random_seed, f, and x.

if nargin < 4
    noise_level = 2e-1;  % The noise level.
end
if nargin < 5
    with_failure = true;  % Whether to fail the function evaluation randomly.
end

% Set the random seed.
orig_rng_state = rng();
rseed = max(0, min(2^32 - 1, random_seed + sum(num2str(f, 16)) + sum(num2str(x, 16), 'all')));
rng(rseed);

% Contaminate f. The value will be further modified below.
f = f * (1 + noise_level * randn);

% Generate a random number to decide how to modify f.
r = 2 * rand - 1;

% Restore the random seed. Do this before the possible invocation of `error`.
rng(orig_rng_state);

% Modify the value of f to make it "tough".
if r > 0.9
    if with_failure
        error('Function evaluation fails!');
    else
        f = NaN;
    end
elseif r > 0.8
    f = NaN;
elseif r > 0.7
    f = Inf;
elseif r < -0.9
    f = -1e30;
end

% `contaminate` ends here
return
