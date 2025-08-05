function [f, g, H] = chrosen(x) % Chained Rosenbrock function

alpha = 4.0;

f = sum((x(1:end-1)-1).^2 + alpha*(x(2:end)-x(1:end-1).^2).^2);  % Function value

n = length(x);

if nargout >= 2
    g=zeros(n,1); % Gradient
    for i=1:n-1
        g(i)   = g(i) + 2*(x(i)-1)+alpha*2*(x(i)^2-x(i+1))*2*x(i);
        g(i+1) = g(i+1) - alpha*2*(x(i)^2-x(i+1));
    end
end

if nargout >= 3
    H=zeros(n,n); % Hessian
    for i=1:n-1
        H(i,i)    =  H(i,i)+2+alpha*2*2*(3*x(i)^2-x(i+1));
        H(i,i+1)  =  H(i,i+1)-alpha*2*2*x(i);
        H(i+1,i)  =  H(i+1,i) -alpha*2*2*x(i);
        H(i+1,i+1)=  H(i+1,i+1)+alpha*2;
    end
end
