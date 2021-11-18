rng(314);
maxdim = 1000;
sigma = 1e6;
nr = 10;
subnr = 10;
for ir = 1 : nr + 3
    ir
    if ir == nr + 1
        m = maxdim;
    else
        m = max(2, floor(rand*maxdim));
    end
    for subir = 1 : subnr + 4
        subir
        if subir == subnr + 1
            n = m;
        else
            n = max(2, floor(rand*m));
        end
        %r = floor(2*rand*min(m, n));
        r = n;
        A = randn(m, r)*diag(1:r)*randn(r, n);
        [Q, R] = qr(A);
        i = ceil(rand*n);
        Rdiag= NaN(min(m,n), 1);
        for i = 1: min(m, n)
            Rdiag(i) = R(i,i);
        end
        [Q, Rdiag] = qrexc(A, Q, Rdiag, i);
        assert(myistriu(Q'*A, 1.0E2*eps), 'Q^T*A is upper triangular');
        assert(norm(diag(Q'*A)-Rdiag)/norm(Rdiag) <= eps, 'Rdiag=diag(R)');
    end
end

function is_triu= myistriu(A, tol)

if (nargin == 2)
    tol_loc = max(tol, tol * max(abs(A)));
else
    tol_loc = ZERO;
end
m = size(A, 1);
n = size(A, 2);
is_triu =true;
for i = 1: min(m, n)
    if (any(abs(A(i + 1:m, i)) > tol_loc))
        is_triu = false;
        break
    end
end
end
