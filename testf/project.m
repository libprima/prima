function y = project(x, V)
if (any(isnan(x)) || any(isnan(V), 'all'))
    y = sum(x) + sum(sum(V));
elseif (any(isinf(V), 'all'))
    V(~isinf(V)) = 0;
    V(isinf(V)) = sign(V(isinf(V)));
    [U,~] = qr(V, 0);
    y = U*(U'*x);
else
    [U,~] = qr(V, 0);
    y = U*(U'*x);
end

tol = 1.0e6*eps;
norm(V'*(x-y))
norm(x-y)
norm(V)
assert(norm(V'*(x-y)) <= max(tol, tol * norm(x - y) * norm(V)), 'X - Y is orthogonal to V');
