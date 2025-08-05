function qdiff = calquad(d, gq, hq, pq, x, xpt, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % CALQUAD calculates
    % QDIFF = Q(X + D) - Q(X)
    % with Q being the quadratic function defined via (GQ, HQ, PQ) by
    % Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
    % where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
    % HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .

    % Local variable

    if debugflag
        funname = 'CALQUAD';

        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(d, n);
        verisize(x, n);
        verisize(gq, n);
        verisize(hq, n, n);
        verisize(pq, npt);
    end

    % The order of calculation seems quite important. The following order seems to work well.
    % 1st order term
    qdiff = d'*gq;
    s = 0.5 * d + x;    % Different from the above version.
    % implicit 2nd-order term
    qdiff = qdiff + sum(pq .* ((s' * xpt)' .* (d' * xpt)'));
    % explicit 2nd-order term
    qdiff = qdiff + s'*(hq * d);
    % The following implementations do not work as well as the above one.
    %qdiff = qdiff + inprod(d, matprod(hq, s))
    %qdiff = qdiff + sum(hq * outprod(s, d))
    %qdiff = qdiff + 0.5_RP*(inprod(d, matprod(hq, s)) + inprod(s, matprod(hq, d)))
end
