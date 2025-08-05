function [idz, bmat, zmat, info] = inith(ij, xpt, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % INITH initializes BMAT and ZMAT.

    bmat = NaN(n, npt + n);
    zmat = NaN(npt, npt - n - 1);

    % Local variables
    funname = 'INITH';

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(ij, 2, npt);
        verisize(bmat, n, npt + n);
        verisize(zmat, npt, npt - n - 1);
    end

    % Set IDZ = 1. It will not be changed in the following.
    idz = 1;

    % Some values to be used for setting BMAT and ZMAT.
    rhobeg = max(abs(xpt(:, 2)));    % Read RHOBEG from XPT.
    rhosq = rhobeg * rhobeg;
    recip = 1 / rhosq;
    reciq = sqrt(0.5) / rhosq;

    % Initialize BMAT and ZMAT to 0.
    bmat = zeros(n, npt + n);
    zmat = zeros(npt, npt - n - 1);

    % Set the nonzero initial elements of BMAT.
    % When NPT >= 2*N + 1, this defines BMAT completely;
    % When NPT <= 2*N, this defines BMAT(1 : NPT - N - 1, :).
    for k = 1:min(npt - n - 1, n)
        bmat(k, k + 1) = 0.5 / rhobeg;
        bmat(k, n + k + 1) = -0.5 / rhobeg;
    end

    % When NPT <= 2*N, set BMAT(NPT - N : N, :)
    for k = (npt - n):n
        bmat(k, 1) = -1 / rhobeg;
        bmat(k, k + 1) = 1 / rhobeg;
        bmat(k, npt + k) = -0.5 * rhosq;
    end

    % Set the nonzero initial elements of ZMAT.
    % When NPT <= 2*N + 1, this defines ZMAT completely;
    % When NPT > 2*N + 1, this defines ZMAT(:, 1 : N).
    for k = 1:min(npt - n - 1, n)
        zmat(1, k) = -reciq - reciq;
        zmat(k + 1, k) = reciq;
        zmat(k + n + 1, k) = reciq;
    end

    % When NPT > 2*N+1, set ZMAT(:, N + 1 : NPT - N - 1).
    for k = (n + 1):(npt - n - 1)
        % I, J, XI, and XJ will be used below.
        i = ij(1, k + n + 1);
        j = ij(2, k + n + 1);
        xi = xpt(i, k + n + 1);
        xj = xpt(j, k + n + 1);
        if xi < 0
            i = i + n;
        end
        if xj < 0
            j = j + n;
        end
        zmat(1, k) = recip;
        zmat(k + n + 1, k) = recip;
        zmat(i + 1, k) = -recip;
        zmat(j + 1, k) = -recip;
    end

    if any(isnan(bmat), 'all') || any(isnan(zmat), 'all')
        info = infos('nan_model');
    else
        info = 0;
    end

end