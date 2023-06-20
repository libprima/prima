function problem = stress_problem(n, problem_type, random_seed)

problem.objective = @(x) chrosen(x);
problem.nonlcon = @(x) nonlcon(x);


function f = chrosen(x)  % the subroutine defining the objective function
alpha = 2;
f = sum((x(1:end-1)-1).^2 + alpha*(x(2:end)-x(1:end-1).^2).^2);
return


function [cineq, ceq] = nonlcon(x)  % the subroutine defining the nonlinear constraints
% The same as fmincon, nonlinear constraints cineq(x) <= 0 and ceq(x) = 0 are specified
% by a function with two returns, the first being cineq and the second being ceq.
cineq = x(2:end) - x(1:end-1).^2;
ceq = x'*x - 1;
return
