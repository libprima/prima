function testeig(n, k)

err = 0;
iter = 0;
for i = 1 : k
td = randn(n,1).*10.^randn(n, 1);
tn = randn(n-1, 1).*10.^randn(n-1, 1);
[it, l1] = eigmin(td, tn);
tridh = spdiags([[tn; 0], td, [0; tn]], -1:1, n, n);
l2 = eigs(tridh, 1, 'smallestreal');

err = max(err, abs(l1-l2)/max(abs(l2), 1));
iter = max(it, iter);
end

err
iter


return


function [iter, eig_min] = eigmin(td, tn)

n = size(td, 1);

maxiter = 50;
tol = 1.0E-8;

pivnew = NaN(n, 1);


piv = -ones(n, 1);
piv(1) = td(1);
for k = 1 : n - 1
    if (piv(k) > 0)
        piv(k + 1) = td(k + 1) - tn(k)^2 / piv(k);
    else
        break
    end
end
if (all(piv >= 0))
        eminub = min(piv);
        eminlb = 0;
else
    eminub = min(td);
    eminlb = -max(abs([0; tn]) + abs(td) + abs([tn; 0]));
end

ksav = 0;
for iter = 1 : maxiter
    if eminub - eminlb <= tol * max(abs(eminlb), abs(eminub))
        break
    end
    eig_min = 0.5 * (eminlb + eminub);

    % The following loop calculates the pivots of the Cholesky factorization of the matrix minus
    % EIG_MIN. All the pivots are positive iff EIG_MIN underestimates the smallest eigenvalue.
    pivnew(1) = td(1) - eig_min;
    for k = 1 : n - 1
        pivnew(k + 1) = td(k + 1) - eig_min - tn(k)^2 / pivnew(k);
    end

    if (all(pivnew > 0))
        piv = pivnew;
        eminlb = eig_min;
        continue
    end

    % We arrive here iff PIVNEW contains nonpositive entries and EIG_MIN is no less than the
    % smallest eigenvalue.
    k = find(~(pivnew > 0), 1, 'first');
    piv(1:k - 1) = pivnew(1:k - 1);
    if (k == ksav && pivksv < 0 && piv(k) - pivnew(k) >= pivnew(k) - pivksv)
            1, iter
            pivksv = 0;
            eminub = (eig_min * piv(k) - eminlb * pivnew(k)) / (piv(k) - pivnew(k));
            %eminub = min(eig_min, (eig_min * piv(k) - eminlb * pivnew(k)) / (piv(k) - pivnew(k)));
    else % K < KSAVE .OR. (K == KSAVE .AND. PIVKSV == 0); note that PIVKSV is always nonpositive.
    %    1
    %    break
        ksav = k;
        pivksv = pivnew(k);  % PIVKSAV <= 0.
        eminub = eig_min;
    end
end

eig_min = eminlb;
