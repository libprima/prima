function [bmat, gq, hq, xbase, xopt, xpt] = shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % SHIFTBASE shifts the base point for XBASE to XBASE + XOPT and updates
    % GQ, HQ, and BMAT accordingly. PQ and ZMAT remain the same after the
    % shifting. See Section 7 of the NEWUOA paper.

    % Local variables
    funname = 'SHIFTBASE';

    % Get and verify the sizes

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(pq, npt);
        verisize(zmat, npt, npt - n - 1);
        verisize(bmat, n, npt + n);
        verisize(gq, n);
        verisize(hq, n, n);
        verisize(xopt, n);
        verisize(xbase, n);
    end

    xoptsq = norm(xopt)^2;
    qxoptq = 0.25 * xoptsq;

    %-------------------------------------------------------------------------%
    %------------------% gq = hessmul(hq, pq, xpt, xopt) + gq %---------------%
    %-----------------------------------% OR %--------------------------------%
    %-% gq = matmul(hq, xopt) + (matmul(xpt, pq * matprod(xopt, xpt)) + gq) %-%
    gq = gq + xpt * (pq .* (xopt' * xpt)');
    gq = gq + hq * xopt;
    %-------------------------------------------------------------------------%

    w1 = (xopt' * xpt)' - 0.5 * xoptsq;
    % W1 equals MATPROD(XPT, XOPT) after XPT is updated as follows.
    xpt = xpt - 0.5 * repmat(xopt, 1, npt);
    %do k = 1, npt
    %    xpt(:, k) = xpt(:, k) - 0.5_RP*xopt
    %end do

    % Update HQ. It has to be done after the above revision to XPT%%%
    xpq = xpt * pq;

    %----------------------------------------------------------------%
    % Implement R2UPDATE properly so that it ensures HQ is symmetric.
    hq = hq + 1 * (xopt * xpq' + xpq * xopt');
    %----------------------------------------------------------------%

    % Make the changes to BMAT that do not depend on ZMAT.
    for k = 1:npt
        bmatk = bmat(:, k);
        w2 = w1(k) * xpt(:, k) + qxoptq * xopt;
        % Implement R2UPDATE properly so that it ensures
        % bmat(:, npt+1:npt+n) is symmetric.
        bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + 1 * (bmatk * w2' + w2 * bmatk');
    end

    % Then the revisions of BMAT that depend on ZMAT are calculated.
    sumz = sum(zmat, 1)';
    for k = 1:(idz - 1)
        %----------------------------------------------------------------------%
        %---------%vlag = qxoptq*sumz(k)*xopt + matprod(xpt, w1*zmat(:, k)) %--%
        vlag = qxoptq * sumz(k) * xopt + xpt * (w1 .* zmat(:, k));
        %----------------------------------------------------------------------%
        bmat(:, 1:npt) = bmat(:, 1:npt) + -1 * (vlag * zmat(:, k)');
        % Implement R1UPDATE properly so that it ensures
        % bmat(:, npt+1:npt+n) is symmetric.
        bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + -1 * (vlag * vlag');
    end
    for k = idz:(npt - n - 1)
        %----------------------------------------------------------------------%
        %---------%vlag = qxoptq*sumz(k)*xopt + matprod(xpt, w1*zmat(:, k)) %--%
        vlag = qxoptq * sumz(k) * xopt + xpt * (w1 .* zmat(:, k));
        %----------------------------------------------------------------------%
        bmat(:, 1:npt) = bmat(:, 1:npt) + 1 * (vlag * zmat(:, k)');
        % Implement R1UPDATE properly so that it ensures
        % bmat(:, npt+1:npt+n) is symmetric.
        bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + 1 * (vlag * vlag');
    end

    %----------------------------------------------------------------%
    % The following instructions complete the shift of XBASE.
    % Recall the we have already subtracted 0.5_RP*XOPT from XPT.
    % Therefore, overall, the new XPT is XPT - XOPT.
    xpt = xpt - 0.5 * repmat(xopt, 1, npt);
    %do k = 1, npt
    %    xpt(:, k) = xpt(:, k) - 0.5_RP*xopt
    %end do

    xbase = xbase + xopt;
    xopt = zeros(n, 1);

end
