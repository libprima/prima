function problem = stress_problem(n, problem_type, random_seed)
%This function generates a test problem for the stress test.
% n: the dimension of the problem
% problem_type:
%   'u' for unconstrained,
%   'b' for bound constrained,
%   'l' for linearly constrained,
%   'n' for nonlinearly constrained
% random_seed: the random seed for the problem generation

% Set the random seed
orig_rng_state = rng();
rng(random_seed);

% Set the starting point
problem.x0 = randn(n, 1);
problem.objective = @chrosen;

% Set the bound constraints
problem.lb = [];
problem.ub = [];
if strcmp(problem_type, 'b') || (ismember(problem_type, {'l', 'n'}) && rand > 0.5)
    problem.lb = problem.x0 - abs(randn(n, 1));
    problem.ub = problem.x0 + abs(randn(n, 1));
end

% Set the linear constraints
problem.Aeq = zeros(0, n);
problem.beq = zeros(0, 1);
problem.Aineq = zeros(0, n);
problem.bineq = zeros(0, 1);
if strcmp(problem_type, 'l') || (strcmp(problem_type, 'n') && rand > 0.5)
    Aeq = hilb(n);
    ind = randperm(n);
    ind = ind(1: floor(n*rand));
    problem.Aeq = Aeq(ind, :);
    problem.beq = problem.Aeq * randn(n, 1);
    e = ones(n, 1);
    Aineq = full(spdiags([e -2*e e], -1:1, n, n));
    ind = randperm(n);
    ind = ind(1: floor(n*rand));
    problem.Aineq = Aineq(ind, :);
    problem.bineq = problem.Aineq * randn(n, 1) + rand(size(problem.Aineq, 1), 1);
end

% Set the nonlinear constraints
problem.nonlcon = [];
if strcmp(problem_type, 'n')
    problem.nonlcon = @nonlcon;
end

% Restore the random seed
rng(orig_rng_state);

% `stress_problem` ends here
return


function f = chrosen(x)  % the subroutine defining the objective function
alpha = 4;
f = sum((x(1:end-1)-1).^2 + alpha*(x(2:end)-x(1:end-1).^2).^2);
return


function [cineq, ceq] = nonlcon(x)  % the subroutine defining the nonlinear constraints
% The same as fmincon, nonlinear constraints cineq(x) <= 0 and ceq(x) = 0 are specified
% by a function with two returns, the first being cineq and the second being ceq.
cineq = x(2:end) - x(1:end-1).^2;
ceq = x'*x - 1;
return
