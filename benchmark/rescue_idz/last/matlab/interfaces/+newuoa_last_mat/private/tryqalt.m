function [itest, gq, hq, pq] = tryqalt(idz, fval, ratio, smat, zmat, itest, gq, hq, pq, debugflag)
    n = length(gq);
    npt = length(pq);
    % TRYQALT tests whether to replace Q by the alternative model, namely the model that minimizes
    % the F-norm of the Hessian subject to the interpolation conditions. It does the replacement
    % when certain criteria are satisfied (i.e., when ITEST = 3). Note that SMAT = BMAT(:, 1:NPT)
    %
    % See Section 8 of the NEWUOA paper.

    % N.B.:
    % GQ, HQ, and PQ should be IN10T(INOUT) instead of IN10T(OUT). According to the Fortran 2018
    % standard, an IN10T(OUT) dummy argument becomes undefined on invocation of the procedure.
    % Therefore, if the procedure does not define such an argument, its value becomes undefined,
    % which is the case for HQ and PQ when ITEST < 3 at exit. In addition, the information in GQ is
    % needed for defining ITEST, so it must be IN10T(INOUT).

    % Local variables
    galt = NaN(length(gq), 1);
    funname = 'TRYQALT';

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(GQ) or SIZE(PQ) is invalid.', funname);
        end
        verisize(fval, npt);
        verisize(smat, n, npt);
        verisize(zmat, npt, npt - n - 1);
        verisize(hq, n, n);
    end

    % In the NEWUOA paper, Powell replaces Q with Q_alt when RATIO <= 0.01 and ||G_alt|| <= 0.1||GQ||
    % hold for 3 consecutive times (eq(8.4)). But Powell's code compares ABS(RATIO) instead of RATIO
    % with 0.01. Here we use RATIO, which is more efficient as observed in in Zhang Zaikun's PhD thesis
    % (Section 3.3.2).
    %if (abs(ratio) > 1.0e-2_RP) then
    if ratio > 1.0e-2
        itest = 0;
    else
        galt = smat * fval;
        if norm(gq)^2 < 1.0e2 * norm(galt)^2
            itest = 0;
        else
            itest = itest + 1;
        end
    end

    % Replace Q with Q_alt when ITEST >= 3.
    if itest >= 3
        gq = galt;
        hq = zeros(n, n);
        fz = (fval' * zmat)';
        fz(1:idz - 1) = -fz(1:idz - 1);
        pq = zmat * fz;
        itest = 0;
    end

end
