maxdim = 1000;
sigma = 1e6;
nr = 10;
subnr = 10;
for ir = 1 : nr + 3
    ir
    if ir == nr + 1
        m = 0;
    elseif ir == nr + 1
        m = 1;
    elseif ir == nr + 3
        m = maxdim;
    else
        m = floor(rand*maxdim);
    end
    for subir = 1 : subnr + 4
        subir
        if subir == subnr + 1
            n = 0;
        elseif subir == subnr + 2
            n = 1;
        elseif subir == subnr + 3
            n = m;
        elseif subnr == subnr + 4
            n = 2*m;
        else
            n = floor(2*rand*m);
        end
        r = floor(2*rand*min(m, n));
        A = sigma*randn(m, r)*diag(sigma.^randn(r, 1))*randn(r, n);
        economy = (rand > 0.5);
        debugging = (rand > 0.2);
        pivote = (rand > 0.5);
        if pivote
            [Q, R, P] = qrfac(A, economy, debugging);
        else
            [Q, R] = qrfac(A, economy, debugging);
        end
    end
end

function [Q, R, P] = qrfac(A, economy, debugging)

pivote = (nargout == 3);
eco = (nargin == 2 && economy);
dbg = (nargin == 3 && debugging);
m = size(A, 1);
n = size(A, 2);

Qfull = eye(m, m);
Rfull = A;
if (pivote)
    P = (1 : n)';
end

for j = 1 : n
    if (pivote)
        [~, k] = max(sum(Rfull(j:m, j:n).^2, 1));
        if (k > 1)
            k = k + j - 1;
            P([j, k]) = P([k, j]);
            Rfull(:, [j, k]) = Rfull(:, [k, j]);
        end
    end
    for i = m : -1: j + 1
        G = planerot(Rfull([j, i], j));
        Rfull([j, i], j) = [hypot(Rfull(j, j), Rfull(i, j)); 0];
        Rfull([j, i], j + 1:n) = G * Rfull([j, i], j + 1:n);
        Qfull(:, [j, i]) = Qfull(:, [j, i]) * G';
    end
end

if (eco && m > n)
    Q = Qfull(:, 1 : n);
    R = Rfull(1 : n, :);
else
    Q = Qfull;
    R = Rfull;
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
