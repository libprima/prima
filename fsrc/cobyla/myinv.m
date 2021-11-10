function B = myinv(A)

n = size(A, 2);

if (istril(A))
    B = zeros(n, n);
    for i = 1 : n
        B(i, i) = 1 / A(i, i);
        B(i, 1:i - 1) = -(A(i, 1:i - 1) / A(i, i))*B(1:i - 1, 1:i - 1);
    end
elseif (istriu(A))
    B = zeros(n, n);
    for i = 1 : n
        B(i, i) = 1 / A(i, i);
        B(1:i - 1, i) = -B(1:i - 1, 1:i - 1)*(A(1:i - 1, i) / A(i, i));
    end
else
    [Q, R, P] = qr(A, 'vector');
    R = R';
    %[Q, R] = qr(A);
    B = zeros(n, n);
    for i = n : -1 : 1
        B(:, i) = (Q(:, i) - B(:, i+1:n)*R(i+1:n, i)) / R(i, i);
    end
    PI(P) = (1:n);
    B = B(:, PI)';
end

tol = 1.0E2*n*eps;
assert(norm(B*A-eye(n)) <= tol* norm(A));
