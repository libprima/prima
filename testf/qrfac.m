function [Q, R, P] = qrfac(A, economy, debugging)

pivote = (nargout == 3);
eco = (nargin == 2 && economy);
dbg = (nargin == 3 && debugging);
m = size(A, 1);
n = size(A, 2);

Qfull = eye(m, m);
T = A';
if (pivote)
    P = (1 : n)';
end

for j = 1 : n
    if (pivote)
        [~, k] = max(sum(T(j:n, j:m).^2, 2));
        if (k > 1)
            k = k + j - 1;
            P([j, k]) = P([k, j]);
            T([j, k], :) = T([k, j], :);
        end
    end
    for i = m : -1: j + 1
        G = planerot(T(j, [j, i])')';
        T(j, [j, i]) = [hypot(T(j, j), T(j, i)); 0];
        T(j + 1:n, [j, i]) = T(j + 1:n, [j, i])*G;
        Qfull(:, [j, i]) = Qfull(:, [j, i]) * G;
    end
end

if (eco && m > n)
    Q = Qfull(:, 1 : n);
    R = T(:, 1 : n)';
else
    Q = Qfull;
    R = T';
end

if (dbg)
    tol = 1.0e1*eps*max(m, n);
    try
        assert(norm(Q'*Q - eye(size(Q, 2))) <= tol, 'The columns of Q are orthogonal');
        assert(istriu(R), 'R is upper triangular');
        if (pivote)
            assert(norm(Q*R - A(:, P)) <= tol * max(1, norm(A)), 'NORM(Q*R - A(:, P)) <= TOL * NORM(A)');
            for j = 1 : min(m, n) - 1
                assert(R(j, j) >= R(j + 1, j + 1), 'R(J, J) >= R(J + 1, J + 1)');
                assert(all(R(j, j)^2 >= sum(R(j : min(m, n), j + 1 : n).^2, 1)), 'R(J, J)^2 >= SUM(R(J : MIN(M, N), J + 1 : N).^2');
            end
        else
            assert(norm(Q*R - A) <= tol * max(1, norm(A)), 'NORM(Q*R - A) <= TOL * NORM(A)');
        end
        if (eco && m > n)
            assert(size(Q, 2) == n && size(R, 1) == n, 'SIZE(Q, 2) == N, SIZE(R, 1) == N');
        end
    catch exception
        keyboard
        if (pivote)
            err = max([norm(Q'*Q - eye(size(Q, 2))), norm(Q*R - A(:, P))/norm(A)])
        else
            err = max([norm(Q'*Q - eye(size(Q, 2))), norm(Q*R - A)/norm(A)])
        end
        rethrow(exception)
    end
end

end

function is_triu = istriu(A)
m = size(A, 1);
n = size(A, 2);
is_triu = true;
for i = 1 : min(m, n)
    if any(abs(A(i + 1 : m, i)) > 0)
        is_triu = false;
        break;
    end
end

end
