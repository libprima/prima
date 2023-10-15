function [f, g] = chrosen(x) % Chained Rosenbrock function
a = randn(500,500); cond(a);
alpha = 4.0;
f = sum((x(2:end)-1).^2 + alpha*(x(2:end).^2 - x(1:end-1)).^2);
if nargout >= 2
    n = length(x);
    g = zeros(n, 1);
    for i = 1:n-1
        g(i+1)   = g(i+1) + 2*(x(i+1)-1)+alpha*2*(x(i+1)^2-x(i))*2*x(i+1);
        g(i) = g(i) - alpha*2*(x(i+1)^2-x(i));
    end
end
