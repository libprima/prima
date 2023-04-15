function d = geostep(idz, knew, kopt, bmat, delbar, xpt, zmat, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % This subroutine finds a step D such that the geometry of the interpolation set is improved when
    % XPT(:, KNEW) is changed to XOPT + D, where XOPT = XPT(:, KOPT)

    % Solver-specific module

    % Local variables
    funname = 'GEOSTEP';

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(bmat, n, npt + n);
        verisize(zmat, npt, npt - n - 1);
    end

    xopt = xpt(:, kopt);    % Read XOPT.

    d = biglag(idz, knew, delbar, bmat, xopt, xpt, zmat, debugflag);

    % ALPHA is the KNEW-th diagonal entry of H.
    zknew = zmat(knew, :)';
    zknew(1:idz - 1) = -zknew(1:idz - 1);
    alpha = zmat(knew, :)*zknew;

    % Calculate VLAG and BETA for D.
    [beta, vlag] = vlagbeta(idz, kopt, bmat, d, xpt, zmat, debugflag);

    % If the cancellation in DENOM is unacceptable, then BIGDEN calculates an alternative model step D.
    if abs(1 + alpha * beta / vlag(knew)^2) <= 0.8
        d = bigden(idz, knew, kopt, bmat, d, xpt, zmat, debugflag);
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = biglag(idz, knew, delbar, bmat, x, xpt, zmat, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % BIGLAG calculates a D by approximately solving
    %
    % max |LFUNC(X + D)|, subject to ||D|| <= DELBAR,
    %
    % where LFUNC is the KNEW-th Lagrange function. See Section 6 of the NEWUOA paper.

    % Local variables
    cf = NaN(5, 1);
    funname = 'BIGLAG';

    % N is the number of variables.
    % NPT is the number of interpolation equations.
    % XPT contains the current interpolation points.
    % BMAT provides the last N ROWs of H.
    % ZMAT and IDZ give a factorization of the first NPT by NPT sub-matrix of H.
    % KNEW is the index of the interpolation point to be dropped.
    % DELBAR is the trust region bound for BIGLAG.
    % D will be set to the step from X to the new point.
    % HCOL, GC, GD, S and W will be used for working space.

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(x, n);
        verisize(bmat, n, npt + n);
        verisize(zmat, npt, npt - n - 1);
    end

    % Set HCOL to the leading NPT elements of the KNEW-th column of H.
    zknew = zmat(knew, :)';
    zknew(1:idz - 1) = -zknew(1:idz - 1);
    hcol = zmat * zknew;

    % Set the unscaled initial direction D. Form the gradient of LFUNC at X, and multiply D by the
    % Hessian of LFUNC.
    d = xpt(:, knew) - x;
    dd = norm(d)^2;

    gd = xpt * (hcol .* (d' * xpt)');

    %----------------------------------------------------------------%
    %-----%gc = bmat(:, knew) + matprod(xpt, hcol*matprod(x, xpt)) %-%
    gc = bmat(:, knew) + xpt * (hcol .* (x' * xpt)');
    %----------------------------------------------------------------%

    % Scale D and GD, with a sign change if required. Set S to another vector in the initial two
    % dimensional subspace.
    gg = norm(gc)^2;
    sp = d'*gc;
    dhd = d'*gd;
    scaling = delbar / sqrt(dd);
    if sp * dhd < 0
        scaling = -scaling;
    end
    t = 0;
    if sp * sp > 0.99 * dd * gg
        t = 1;
    end
    tau = scaling * (abs(sp) + 0.5 * scaling * abs(dhd));
    if gg * (delbar * delbar) < 1.0e-2 * tau * tau
        t = 1;
    end
    d = scaling * d;
    gd = scaling * gd;
    s = gc + t * gd;

    % Begin the iteration by overwriting S with a vector that has the required length and direction,
    % except that termination occurs if the given D and S are nearly parallel.
    for iterc = 1:n
        dd = norm(d)^2;
        sp = d'*s;
        ss = norm(s)^2;
        if dd * ss - sp * sp <= 1.0e-8 * dd * ss
            break;
        end
        denom = sqrt(dd * ss - sp * sp);
        s = (dd * s - sp * d) / denom;

        w = xpt * (hcol .* (s' * xpt)');

        % Calculate the coefficients of the objective function on the circle, beginning with the
        % multiplication of S by the second derivative matrix.
        cf(1) = s'*w;
        cf(2) = d'*gc;
        cf(3) = s'*gc;
        cf(4) = d'*gd;
        cf(5) = s'*gd;
        cf(1) = 0.5 * cf(1);
        cf(4) = 0.5 * cf(4) - cf(1);

        % Seek the value of the angle that maximizes |TAU|.
        taubeg = cf(1) + cf(2) + cf(4);
        taumax = taubeg;
        tauold = taubeg;
        isav = 0;
        iu = 49;
        unitang = (2 * pi) / (iu + 1);

        for i = 1:iu
            angle = i * unitang;
            cth = cos(angle);
            sth = sin(angle);
            tau = cf(1) + (cf(2) + cf(4) * cth) * cth + (cf(3) + cf(5) * cth) * sth;
            if abs(tau) > abs(taumax)
                taumax = tau;
                isav = i;
                taua = tauold;
            elseif i == isav + 1
                taub = tau;
            end
            tauold = tau;
        end

        if isav == 0
            taua = tau;
        end
        if isav == iu
            taub = taubeg;
        end
        if abs(taua - taub) > 0
            taua = taua - taumax;
            taub = taub - taumax;
            step = 0.5 * (taua - taub) / (taua + taub);
        else
            step = 0;
        end
        angle = unitang * (isav + step);

        % Calculate the new D and GD. Then test for convergence.
        cth = cos(angle);
        sth = sin(angle);
        tau = cf(1) + (cf(2) + cf(4) * cth) * cth + (cf(3) + cf(5) * cth) * sth;

        dold = d;
        d = cth * d + sth * s;
        % Exit in case of Inf/NaN in D.
        if isinf(sum(abs(d)))
            d = dold;
            break;
        end

        gd = cth * gd + sth * w;
        s = gc + gd;
        if abs(tau) <= 1.1 * abs(taubeg)
            break;
        end
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = bigden(idz, knew, kopt, bmat, d0, xpt, zmat, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % BIGDEN calculates a D by approximately solving
    %
    % max |SIGMA(XOPT + D)|, subject to ||D|| <= DELBAR,
    %
    % where SIGMA is the denominator sigma in the updating formula (4.11)--(4.12) for H, which is the
    % inverse of the coefficient matrix for the interplolation system (see (3.12)). Indeed, each column
    % of H corresponds to a Lagrange basis function of the interpolation problem.  See Section 6 of the
    % NEWUOA paper.
    % N.B.:
    % In Powell's code, BIGDEN calculates also the VLAG and BETA for the selected D. Here, to reduce the
    % coupling of code, we return only D but computes VLAG and BETA outside by calling VLAGBETA. This
    % does not change the mathematics, but the computed VLAG (BETA) will be numerically different due to
    % rounding errors.

    % Local variable
    den = NaN(9, 1);
    denex = NaN(9, 1);
    par = NaN(9, 1);
    prod = NaN(size(xpt, 2) + size(xpt, 1), 5);
    w = NaN(size(xpt, 2) + size(xpt, 1), 5);
    funname = 'BIGDEN';

    % N is the number of variables.
    % NPT is the number of interpolation equations.
    % X is the best interpolation point so far.
    % XPT contains the current interpolation points.
    % BMAT provides the last N ROWs of H.
    % ZMAT and IDZ give a factorization of the first NPT by NPT sub-matrix of H.
    % NDIM is the second dimension of BMAT and has the value NPT + N.
    % KOPT is the index of the optimal interpolation point.
    % KNEW is the index of the interpolation point to be dropped.
    % D will be set to the step from X to the new point, and on entry it should be the D that was
    % calculated by the last call of BIGLAG. The length of the initial D provides a trust region bound
    % on the final D.
    %
    % D is calculated in a way that should provide a denominator with a large modulus in the updating
    % formula when the KNEW-th interpolation point is shifted to the new position X + D.

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(bmat, n, npt + n);
        verisize(d0, n);
        verisize(zmat, npt, npt - n - 1);
    end

    x = xpt(:, kopt);    % For simplicity, use x to denote XOPT.

    % Store the first NPT elements of the KNEW-th column of H in HCOL.
    zknew = zmat(knew, :)';
    zknew(1:idz - 1) = -zknew(1:idz - 1);
    hcol = zmat * zknew;
    alpha = hcol(knew);

    % The initial search direction D is taken from the last call of BIGLAG, and the initial S is set
    % below, usually to the direction from X to X_KNEW, but a different direction to an interpolation
    % point may be chosen, in order to prevent S from being nearly parallel to D.
    d = d0;
    dd = norm(d)^2;
    s = xpt(:, knew) - x;
    ds = d'*s;
    ss = norm(s)^2;
    xsq = norm(x)^2;

    if ~(ds * ds <= 0.99 * dd * ss)
        % `.NOT. (A <= B)` differs from `A > B`.  The former holds iff A > B or {A, B} contains NaN.
        dtest = ds * ds / ss;
        xptemp = xpt - repmat(x, 1, npt);
        %----------------------------------------------------------------%
        %---------%dstemp = matprod(d, xpt) - inprod(x, d) %-------------%
        dstemp = (d' * xptemp)';
        %----------------------------------------------------------------%
        sstemp = sum((xptemp).^2, 1)';

        dstemp(kopt) = 2 * ds + 1;
        sstemp(kopt) = ss;
        [~, k] = min(dstemp .* dstemp ./ sstemp);
        % K can be 0 due to NaN. In that case, set K = KNEW. Otherwise, memory errors will occur.
        if k == 0
            k = knew;
        end
        if (~(dstemp(k) * dstemp(k) / sstemp(k) >= dtest)) && k ~= kopt
            % `.NOT. (A >= B)` differs from `A < B`.  The former holds iff A < B or {A, B} contains NaN.
            % Although unlikely, if NaN occurs, it may happen that K = KOPT.
            s = xpt(:, k) - x;
            ds = dstemp(k);
            ss = sstemp(k);
        end
    end

    ssden = dd * ss - ds * ds;
    densav = 0;

    % Begin the iteration by overwriting S with a vector that has the required length and direction.
    for iterc = 1:n
        s = (1 / sqrt(ssden)) * (dd * s - ds * d);
        xd = x'*d;
        xs = x'*s;

        % Set the coefficients of the first two terms of BETA.
        tempa = 0.5 * xd * xd;
        tempb = 0.5 * xs * xs;
        den(1) = dd * (xsq + 0.5 * dd) + tempa + tempb;
        den(2) = 2 * xd * dd;
        den(3) = 2 * xs * dd;
        den(4) = tempa - tempb;
        den(5) = xd * xs;
        den(6:9) = 0;

        % Put the coefficients of WCHECK in W.
        for k = 1:npt
            tempa = xpt(:, k)'*d;
            tempb = xpt(:, k)'*s;
            tempc = xpt(:, k)'*x;
            w(k, 1) = 0.25 * (tempa * tempa + tempb * tempb);
            w(k, 2) = tempa * tempc;
            w(k, 3) = tempb * tempc;
            w(k, 4) = 0.25 * (tempa * tempa - tempb * tempb);
            w(k, 5) = 0.5 * tempa * tempb;
        end
        w(npt + 1:npt + n, 1:5) = 0;
        w(npt + 1:npt + n, 2) = d;
        w(npt + 1:npt + n, 3) = s;

        % Put the coefficents of THETA*WCHECK in PROD.
        for jc = 1:5
            wz = (w(1:npt, jc)' * zmat)';
            wz(1:idz - 1) = -wz(1:idz - 1);
            prod(1:npt, jc) = zmat * wz;

            nw = npt;
            if jc == 2 || jc == 3
                prod(1:npt, jc) = prod(1:npt, jc) + (w(npt + 1:npt + n, jc)' * bmat(:, 1:npt))';
                nw = npt + n;
            end
            prod(npt + 1:npt + n, jc) = bmat(:, 1:nw) * w(1:nw, jc);
        end

        % Include in DEN the part of BETA that depends on THETA.
        for k = 1:(npt + n)
            par(1:5) = 0.5 * prod(k, 1:5)' .* w(k, 1:5)';
            den(1) = den(1) - par(1) - sum(par(1:5));
            tempa = prod(k, 1) * w(k, 2) + prod(k, 2) * w(k, 1);
            tempb = prod(k, 2) * w(k, 4) + prod(k, 4) * w(k, 2);
            tempc = prod(k, 3) * w(k, 5) + prod(k, 5) * w(k, 3);
            den(2) = den(2) - tempa - 0.5 * (tempb + tempc);
            den(6) = den(6) - 0.5 * (tempb - tempc);
            tempa = prod(k, 1) * w(k, 3) + prod(k, 3) * w(k, 1);
            tempb = prod(k, 2) * w(k, 5) + prod(k, 5) * w(k, 2);
            tempc = prod(k, 3) * w(k, 4) + prod(k, 4) * w(k, 3);
            den(3) = den(3) - tempa - 0.5 * (tempb - tempc);
            den(7) = den(7) - 0.5 * (tempb + tempc);
            tempa = prod(k, 1) * w(k, 4) + prod(k, 4) * w(k, 1);
            den(4) = den(4) - tempa - par(2) + par(3);
            tempa = prod(k, 1) * w(k, 5) + prod(k, 5) * w(k, 1);
            tempb = prod(k, 2) * w(k, 3) + prod(k, 3) * w(k, 2);
            den(5) = den(5) - tempa - 0.5 * tempb;
            den(8) = den(8) - par(4) + par(5);
            tempa = prod(k, 4) * w(k, 5) + prod(k, 5) * w(k, 4);
            den(9) = den(9) - 0.5 * tempa;
        end

        par(1:5) = 0.5 * prod(knew, 1:5)'.^2;
        denex(1) = alpha * den(1) + par(1) + sum(par(1:5));
        tempa = 2 * prod(knew, 1) * prod(knew, 2);
        tempb = prod(knew, 2) * prod(knew, 4);
        tempc = prod(knew, 3) * prod(knew, 5);
        denex(2) = alpha * den(2) + tempa + tempb + tempc;
        denex(6) = alpha * den(6) + tempb - tempc;
        tempa = 2 * prod(knew, 1) * prod(knew, 3);
        tempb = prod(knew, 2) * prod(knew, 5);
        tempc = prod(knew, 3) * prod(knew, 4);
        denex(3) = alpha * den(3) + tempa + tempb - tempc;
        denex(7) = alpha * den(7) + tempb + tempc;
        tempa = 2 * prod(knew, 1) * prod(knew, 4);
        denex(4) = alpha * den(4) + tempa + par(2) - par(3);
        tempa = 2 * prod(knew, 1) * prod(knew, 5);
        denex(5) = alpha * den(5) + tempa + prod(knew, 2) * prod(knew, 3);
        denex(8) = alpha * den(8) + par(4) - par(5);
        denex(9) = alpha * den(9) + prod(knew, 4) * prod(knew, 5);

        % Seek the value of the angle that maximizes the |DENOM|.
        denom = denex(1) + denex(2) + denex(4) + denex(6) + denex(8);
        denold = denom;
        denmax = denom;
        isav = 0;
        iu = 49;
        unitang = (2 * pi) / (iu + 1);
        par(1) = 1;
        for i = 1:iu
            angle = i * unitang;
            par(2) = cos(angle);
            par(3) = sin(angle);
            for j = 4:2:8
                par(j) = par(2) * par(j - 2) - par(3) * par(j - 1);
                par(j + 1) = par(2) * par(j - 1) + par(3) * par(j - 2);
            end
            denomold = denom;
            denom = denex(1:9)'*par(1:9);
            if abs(denom) > abs(denmax)
                denmax = denom;
                isav = i;
                dena = denomold;
            elseif i == isav + 1
                denb = denom;
            end
        end
        if isav == 0
            dena = denom;
        end
        if isav == iu
            denb = denold;
        end
        if abs(dena - denb) > 0
            dena = dena - denmax;
            denb = denb - denmax;
            step = 0.5 * (dena - denb) / (dena + denb);
        else
            step = 0;
        end
        angle = unitang * (isav + step);

        % Calculate the new parameters of the denominator, the new VLAG vector and the new D. Then test
        % for convergence.
        par(2) = cos(angle);
        par(3) = sin(angle);
        for j = 4:2:8
            par(j) = par(2) * par(j - 2) - par(3) * par(j - 1);
            par(j + 1) = par(2) * par(j - 1) + par(3) * par(j - 2);
        end

        %beta = inprod(den(1:9), par(1:9))  % Not needed since we do not return BETA.

        denmax = denex(1:9)'*par(1:9);

        vlag = prod(:, 1:5) * par(1:5);

        tau = vlag(knew);

        dold = d;
        d = par(2) * d + par(3) * s;
        % Exit in case of Inf/NaN in D.
        if isinf(sum(abs(d)))
            d = dold;
            break;
        end

        dd = norm(d)^2;
        xnew = x + d;
        dxn = d'*xnew;
        xnsq = norm(xnew)^2;

        if iterc > 1
            densav = max(densav, denold);
        end
        if abs(denmax) <= 1.1 * abs(densav)
            break;
        end
        densav = denmax;

        % Set S to 0.5_RP the gradient of the denominator with respect to D.
        s = tau * bmat(:, knew) + alpha * (dxn * x + xnsq * d - vlag(npt + 1:npt + n));
        v = (xnew' * xpt)';
        v = (tau * hcol - alpha * vlag(1:npt)) .* v;
        %------------------------------------------------------%
        %---------%s = s + matprod(xpt, v) %-------------------%
        s = s + xpt * v;
        %------------------------------------------------------%

        ss = norm(s)^2;
        ds = d'*s;
        ssden = dd * ss - ds * ds;
        if ssden < 1.0e-8 * dd * ss
            break;
        end
    end

    %vlag(kopt) = vlag(kopt) + 1  % Note needed since we do not return VLAG.

end
