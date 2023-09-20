function f = chrosen(x) % Chained Rosenbrock function
f = sum((x(1:end-1)-1).^2 + 4*(x(2:end)-x(1:end-1).^2).^2);
end
