
clc
clear all;
close all;

problem = struct();
problem.objective = @ChenxianyunFunc;
problem.x0 = [10,10];
problem.lb = [0,0];
problem.ub = [50,50];
problem.options.maxfun = 1000000;
problem.options.rhoend = 0;
[a, b, c, d] = cobyla(problem)


function z = ChenxianyunFunc(x)
% 答案：
% 目标函数值(最小): -715.775061146149
% x1: 20.4786412351456
% x2: 25.066028561025
x1=x(1);
x2=x(2);
z = -40*x1*cos(x2)-35*x2*sin(x1)+15*x1^2-sin(x1)*20+12*x2^2-25*x1*x2*sin(x1*x2^2*cos(x1+x2));
end
