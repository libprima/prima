function [cineq, ceq] = cfail(x)
cineq = [log(x), -1];
ceq = sqrt(x);
if abs(x) < 0.5
    error('C fail at %f', x);
end
