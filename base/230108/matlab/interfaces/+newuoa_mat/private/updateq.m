function [gq, hq, pq] = updateq(idz, knew, bmatknew, fqdiff, zmat, xptknew, gq, hq, pq, debugflag)
    n = length(gq);
    npt = length(pq);
    % UPDATEQ updates GQ, HQ, and PQ when XPT(:, KNEW) is replaced by XNEW.
    % See Section 4 of the NEWUOA paper.
    %
    % FQDIFF = [F(XNEW) - F(XOPT)] - [Q(XNEW) - Q(XOPT)] = MODERR

    % Local variables
    funname = 'UPDATEQ';

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(GQ) or SIZE(PQ) is invalid.', funname);
        end
        verisize(bmatknew, n);
        verisize(zmat, npt, npt - n - 1);
        verisize(xptknew, n);
        verisize(hq, n, n);
    end

    %----------------------------------------------------------------%
    % Implement R1UPDATE properly so that it ensures HQ is symmetric.
    hq = hq + pq(knew) * (xptknew * xptknew');
    %----------------------------------------------------------------%

    % Update the implicit part of second derivatives.
    fqdz = fqdiff * zmat(knew, :)';
    fqdz(1:idz - 1) = -fqdz(1:idz - 1);
    pq(knew) = 0;
    %----------------------------------------------------------------%
    %pq = pq + matprod(zmat, fqdz) %---------------------------------%
    pq = pq + zmat * fqdz;
    %----------------------------------------------------------------%

    % Update the gradient.
    gq = gq + fqdiff * bmatknew;

end
