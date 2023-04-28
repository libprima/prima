function knew = setdrop_tr(idz, kopt, beta, delta, ratio, rho, vlag, xopt, xpt, zmat, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % This subroutine sets KNEW to the index of the interpolation point to be deleted AFTER A TRUST
    % REGION STEP. KNEW will be set in a way ensuring that the geometry of XPT is "optimal" after
    % XPT(:, KNEW) is replaced by XNEW = XOPT + D, where D is the trust-region step.  Note that the
    % information of XNEW is included in VLAG and BETA, which are calculated according to D.
    %
    % N.B.: At the entry of this function is invoked, XOPT may differ from XPT(:, KOPT), because XOPT is
    % updated but KOPT is not. See NEWUOB for details.

    % Local variables
    funname = 'SETDROP_TR';

    % Get and verify the sizes

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(vlag, npt);
        verisize(xopt, n);
        verisize(zmat, npt, npt - n - 1);
    end

    rhosq = max(0.1 * delta, rho)^2;
    hdiag = -sum(zmat(:, 1:idz - 1).^2, 2) + sum(zmat(:, idz:npt - n - 1).^2, 2);
    xdsq = sum((xpt - repmat(xopt, 1, npt)).^2, 1)';
    sigma = abs(beta * hdiag + vlag(1:npt).^2);
    sigma = sigma .* max(xdsq / rhosq, 1).^3;
    if ratio <= 0
        % If the new F is not better than the current FVAL(KOPT), we set SIGMA(KOPT) = -1 to prevent
        % KNEW from being KOPT.
        sigma(kopt) = -1;
    end
    if ratio > 0 || (max(sigma) > 1 && ~any(isnan(sigma)))
        % N.B.:
        % 1. In Powell's code, the above condition is (RATIO > 0 .OR. MAXVAL(SIGMA) > 1).
        % 2. THEORETICALLY speaking, with Powell's condition, when RATIO > 0 (i.e., the new F is less
        % than the current FVAL(KOPT)), the following line sets KNEW > 0, ensuring XNEW to be included
        % in XPT. However, KNEW may turn out 0 due to NaN in SIGMA: for GFortran, KNEW = 0 if SIGMA
        % contains only NaN, yet other compilers/languages may behave differently.
        % 3. With the new condition, when SIGMA contains NaN (can happen in bad-conditioned problems,
        % although rarely), we explicitly set KNEW = 0. Consequently, NEWUOA will check whether to take
        % a geometry step or reduce RHO.
        [~, knew] = max(sigma);
    else
        knew = 0;
    end
    % It is tempting to take the function value into consideration when defining KNEW, for example,
    % set KNEW so that FVAL(KNEW) = MAX(FVAL) as long as F(XNEW) < MAX(FVAL), unless there is a better
    % choice. However, this is not a good idea, because the definition of KNEW should benefit the
    % quality of the model that interpolates f at XPT. A set of points with low function values is not
    % necessarily a good interpolation set. In contrast, a good interpolation set needs to include
    % points with relatively high function values; otherwise, the interpolant will unlikely reflect the
    % landscape of the function sufficiently.
end
