%ROSENBROCK_EXAMPLE illustrates how to use pdfo.
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

fprintf('\nMinimize the chained Rosenbrock function with three variables subject to various constraints:\n');
x0 = [0; 0; 0];  % starting point

fprintf('\n1. Nonlinear constraints --- ||x||_2^2 = 1, x(i)^2 >= x(i+1) >= 0.5*x(i) >= 0 for i = 1, 2:\n');
% linear inequality constraints A*x <= b
A = [0.5, -1, 0; 0, 0.5, -1];
b = [0; 0];
% linear equality constraints Aeq*x = beq
Aeq = [];
beq = [];
% bound constraints lb <= x <= ub
lb = [0; 0; 0];
ub = [];  % ub = [inf; inf; inf] works equally well
% nonlinear constraints
nonlcon = @nlc;  % see function nlc given below
% The following syntax is identical to fmincon:
[x, fx, exitflag, output] = pdfo(@chrosen, x0, A, b, Aeq, beq, lb, ub, nonlcon)
% Alternatively, the problem can be passed to pdfo as a structure:
%p.objective = @chrosen; p.x0 = x0; p.Aineq = A; p.bineq = b; p.lb = lb; p.nonlcon = @nlc;
%[x, fx, exitflag, output] = pdfo(p)

fprintf('\n2. Linear constraints --- sum(x) = 1, x(i+1) <= x(i) <= 1 for i = 1, 2:\n');
A = [-1, 1, 0; 0, -1, 1];
b = [0; 0];
Aeq = [1, 1, 1];
beq = 1;
ub = [1; 1; 1];
[x, fx, exitflag, output] = pdfo(@chrosen, x0, A, b, Aeq, beq, [], ub)

fprintf('\n3. Bound constraints --- -0.5 <= x(1) <= 0.5, 0 <= x(2) <= 0.25:\n');
lb = [-0.5; 0; -inf];
ub = [0.5; 0.25; inf];
[x, fx, exitflag, output] = pdfo(@chrosen, x0, [], [], [], [], lb, ub)

fprintf('\n4. No constraints:\n');
[x, fx, exitflag, output] = pdfo(@chrosen, x0)


function f = chrosen(x)  % the subroutine defining the objective function
f = sum((x(1:end-1)-1).^2 + 4*(x(2:end)-x(1:end-1).^2).^2);
end


function [cineq, ceq] = nlc(x)  % the subroutine defining the nonlinear constraints
% The same as fmincon, nonlinear constraints cineq(x) <= 0 and ceq(x) = 0 are specified
% by a function with two returns, the first being cineq and the second being ceq.
cineq = x(2:end) - x(1:end-1).^2;
ceq = x'*x - 1;
end
