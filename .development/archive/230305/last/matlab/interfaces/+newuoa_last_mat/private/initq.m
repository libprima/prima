function [gq, hq, pq, info] = initq(ij, fval, xpt, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % INITQ initializes the quadratic model, which is represented by
    % (GQ, HQ, PQ) in the following way:
    % the gradient of the model at XBASE is GQ;
    % the Hessian of the model is
    % HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.

    gq = NaN(n, 1);
    hq = NaN(n, n);
    pq = NaN(npt, 1);

    % Local variables
    funname = 'INITQ';

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(fval, npt);
        verisize(ij, 2, npt);
        verisize(gq, n);
        verisize(hq, n, n);
        verisize(pq, npt);
    end

    gq = zeros(n, 1);
    hq = zeros(n, n);
    pq = zeros(npt, 1);    % We will not update PQ. It is 0 at return.

    rhobeg = max(abs(xpt(:, 2)));    % Read RHOBEG from XPT.
    fbeg = fval(1);

    % Set GQ by forward difference.
    gq(1:n) = (fval(2:n + 1) - fbeg) / rhobeg;
    % If possible, revise GQ to central difference.
    k = min(npt - n - 1, n);
    gq(1:k) = 0.5 * (gq(1:k) + (fbeg - fval(n + 2:n + 1 + k)) / rhobeg);

    % Set the diagonal of HQ by 2nd-order central finite difference.
    for k = 1:min(npt - n - 1, n)
        hq(k, k) = ((fval(k + 1) - fbeg) / rhobeg - (fbeg - fval(k + n + 1)) / rhobeg) / rhobeg;
    end
    % When NPT > 2*N + 1, set the off-diagonal entries of HQ.
    for k = (2 * n + 2):npt
        % I, J, XI, and XJ will be used below.
        i = ij(1, k);
        j = ij(2, k);
        xi = xpt(i, k);
        xj = xpt(j, k);
        if xi * xpt(i, i + 1) > 0
            fi = fval(i + 1);
        else
            fi = fval(i + n + 1);
        end
        if xj * xpt(j, j + 1) > 0
            fj = fval(j + 1);
        else
            fj = fval(j + n + 1);
        end
        % With the XI, XJ, FI, and FJ found above, we have
        % FVAL(K) = F(XBASE + XI + XJ),
        % FI = F(XBASE + XI),
        % FJ = F(XBASE + XJ).
        % Thus the HQ(I, J) defined below approximates
        % frac{partial^2}{partial X_I partial X_J} F(XBASE)
        hq(i, j) = (fbeg - fi - fj + fval(k)) / (xi * xj);
        hq(j, i) = hq(i, j);
    end

    if any(isnan(gq)) || any(isnan(hq), 'all') || any(isnan(pq))
        info = infos('nan_model');
    else
        info = 0;
    end

end
