function f = chrosen(x) % Chained Rosenbrock function
alpha = 4.0;
f = sum((x(1:end-1)-1).^2 + alpha*(x(2:end)-x(1:end-1).^2).^2);
end
