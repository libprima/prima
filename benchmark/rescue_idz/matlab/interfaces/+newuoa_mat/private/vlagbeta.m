function [beta, vlag] = vlagbeta(idz, kopt, bmat, d, xpt, zmat, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % VLAGBETA is calculates VLAG = Hw and BETA for a given step D.
    % See (4.11)--(4.12) of the NEWUOA paper.

    vlag = NaN(npt + n, 1);

    % Local variables
    funname = 'VLAGBETA';

    % Get and verify the sizes

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(bmat, n, npt + n);
        verisize(zmat, npt, npt - n - 1);
        verisize(d, n);
        verisize(vlag, n + npt);
    end

    xopt = xpt(:, kopt);    % Read XOPT.

    %----------------------------------------------------------------------%
    % This is the one of the two places where WCHECK is calculated,
    % the other one being BIGDEN (now removed).
    % WCHECK contains the first NPT entries of (w-v) for the vectors
    % w and v defined in eq(4.10) and eq(4.24) of the NEWUOA paper,
    % and also hat{w} in eq(6.5) of
    %
    % M. J. D. Powell, Least Frobenius norm updating of quadratic
    % models that satisfy interpolation conditions. Math. Program.,
    % 100:183--215, 2004
    wcheck = (d' * xpt)';
    wcheck = wcheck .* (0.5 * wcheck + (xopt' * xpt)');
    %----------------------------------------------------------------------%

    vlag(1:npt) = (d' * bmat(:, 1:npt))';

    wz = (wcheck' * zmat)';
    wzsav = wz;
    wz(1:idz - 1) = -wz(1:idz - 1);
    beta = -(wzsav'*wz);
    %----------------------------------------------------------------------%
    %-----%vlag(1 : npt) = vlag(1 : npt) + matprod(zmat, wz) %-------------%
    vlag(1:npt) = vlag(1:npt) + zmat * wz;
    %----------------------------------------------------------------------%

    vlag(kopt) = vlag(kopt) + 1;    % The calculation of VLAG(1:NPT) finishes.

    bw = bmat(:, 1:npt) * wcheck;
    %----------------------------------------------------------------------%
    %vlag(npt + 1 : npt + n) = bw + matprod(d, bmat(:, npt + 1 : npt + n)) %
    vlag(npt + 1:npt + n) = bw + (d' * bmat(:, npt + 1:npt + n))';
    % The calculation of VLAG finishes.
    %----------------------------------------------------------------------%

    %----------------------------------------------------------------------%
    %-----%bwvd = inprod(bw + vlag(npt + 1 : npt + n), d) %----------------%
    bwvd = dot(bw + vlag(npt + 1:npt + n), d);
    %----------------------------------------------------------------------%

    dx = d'*xopt;
    dsq = norm(d)^2;
    xoptsq = norm(xopt)^2;

    % The final value of BETA is calculated as follows.
    beta = dx * dx + dsq * (xoptsq + dx + dx + 0.5 * dsq) + beta - bwvd;

end
