function [crvmin, s, info] = trsapp(delta, gq, hq, pq, tol, x, xpt, debugflag)
    n = size(xpt, 1);
    npt = size(xpt, 2);
    % TRSAPP finds an approximate solution to the N-dimensional trust
    % region subproblem
    %
    % min <X+S, GQ> + 0.5*<X+S, HESSIAN*(X+S)> s.t. ||S|| <= DELTA
    %
    % Note that the HESSIAN here is the sum of an explicit part HQ and
    % an implicit part (PQ, XPT):
    %
    % HESSIAN = HQ + sum_K=1^NPT PQ(K)*XPT(:, K)*XPT(:, K)' .
    %
    % At return, S will be the approximate solution. CRVMIN will be
    % set to the least curvature of HESSIAN along the conjugate
    % directions that occur, except that it is set to 0 if S goes
    % all the way to the trust region boundary. QRED is the reduction
    % of Q achieved by S. INFO is an exit flag:
    % INFO = 0: an approximate solution satisfying one of the
    % following conditions is found:
    % 1. ||G+HS||/||GBEG|| <= TOL,
    % 2. ||S|| = DELTA and <S, -(G+HS)> >= (1 - TOL)*||S||*||G+HS||,
    % where TOL is a tolerance that is set to 1e-2 in NEWUOA.
    % INFO = 1: the iteration is reducing Q only slightly;
    % INFO = 2: the maximal number of iterations is attained;
    % INFO = -1: too much rounding error to continue

    % The calculation of S begins with the truncated conjugate
    % gradient method. If the boundary of the trust region is reached,
    % then further changes to S may be made, each one being in the 2D
    % space spanned by the current S and the corresponding gradient of
    % Q. Thus S should provide a substantial reduction to Q within the
    % trust region.
    %
    % See Section 5 of the NEWUOA paper.

    s = NaN(n, 1);

    % Local variables
    funname = 'TRSAPP';

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(XPT) is invalid.', funname);
        end
        verisize(x, n);
        verisize(gq, n);
        verisize(hq, n, n);
        verisize(pq, npt);
        verisize(s, n);
    end

    s = zeros(n, 1);
    crvmin = 0;
    qred = 0;
    info = 2;    % Default exit flag is 2, i.e., itermax is attained

    % Prepare for the first line search.
    %----------------------------------------------------------------%
    %-----%hx = hessmul(hq, pq, xpt, x) %----------------------------%
    hx = xpt * (pq .* (x' * xpt)') + hq * x;
    %----------------------------------------------------------------%
    g = gq + hx;
    gg = norm(g)^2;
    gg0 = gg;
    d = -g;
    dd = gg;
    ds = 0;
    ss = 0;
    hs = zeros(length(x), 1);
    delsq = delta * delta;
    itermax = n;

    twod_search = false;

    % The truncated-CG iterations.
    %
    % The iteration will be terminated in 4 possible cases:
    % 1. the maximal number of iterations is attained;
    % 2. QADD <= TOL*QRED or ||G|| <= TOL*||GBEG||, where QADD is the
    %    reduction of Q due to the latest CG step, QRED is the
    %    reduction of Q since the beginning until the latest CG step,
    %    G is the current gradient, and GBEG is the initial gradient;
    %    see (5.13) of the NEWUOA paper;
    % 3. DS <= 0
    % 4. ||S|| = DELTA, i.e., CG path cuts the trust region boundary.
    %
    % In the 4th case, twod_search will be set to true, meaning that S
    % will be improved by a sequence of two-dimensional search, the
    % two-dimensional subspace at each iteration being span(S, -G).
    for iter = 1:itermax
        % Check whether to exit due to small GG
        if gg <= (tol^2) * gg0
            info = 0;
            break;
        end
        % Set BSTEP to the step length such that ||S+BSTEP*D|| = DELTA
        bstep = (delsq - ss) / (ds + sqrt(ds * ds + dd * (delsq - ss)));
        %----------------------------------------------------------------%
        %-----%hd = hessmul(hq, pq, xpt, d) %----------------------------%
        hd = xpt * (pq .* (d' * xpt)') + hq * d;
        %----------------------------------------------------------------%
        dhd = d'*hd;

        % Set the step-length ALPHA and update CRVMIN.
        if dhd <= 0
            alpha = bstep;
        else
            alpha = min(bstep, gg / dhd);
            if iter == 1
                crvmin = dhd / dd;
            else
                crvmin = min(crvmin, dhd / dd);
            end
        end
        % QADD is the reduction of Q due to the new CG step.
        qadd = alpha * (gg - 0.5 * alpha * dhd);
        % QRED is the reduction of Q up to now.
        qred = qred + qadd;
        % QADD and QRED will be used in the 2D minimization if any.

        % Update S, HS, and GG.
        sold = s;
        s = s + alpha * d;
        ss = norm(s)^2;
        hs = hs + alpha * hd;
        ggsav = gg;        % Gradient norm square before this iteration
        gg = norm(g + hs)^2;        % Current gradient norm square
        % We may record g+hs for later usage:
        % gnew = g + hs
        % Note that we should NOT set g = g + hs, because g contains
        % the gradient of Q at x.

        % Check whether to exit. This should be done after updating HS
        % and GG, which will be used for the 2D minimization if any.
        if alpha >= bstep || ss >= delsq
            % CG path cuts the boundary. Set CRVMIN to 0.
            crvmin = 0;
            % The only possibility that twod_search is true.
            twod_search = true;
            break;
        end

        % Check whether to exit due to small QADD
        if qadd <= tol * qred
            info = 1;
            break;
        end

        % Exit in case of Inf/NaN in S.
        if isinf(sum(abs(s)))
            s = sold;
            info = -1;
            break;
        end

        % Prepare for the next CG iteration.
        d = (gg / ggsav) * d - g - hs;        % CG direction
        dd = norm(d)^2;
        ds = d'*s;
        if ds <= 0
            % DS is positive in theory.
            info = -1;
            break;
        end
    end

    if ss <= 0 || isnan(ss)
        % This may occur for ill-conditioned problems due to rounding.
        info = -1;
        twod_search = false;
    end

    if twod_search
        % At least 1 iteration of 2D minimization
        itermax = max(1, itermax - iter);
    else
        itermax = 0;
    end

    % The 2D minimization
    for iter = 1:itermax
        if gg <= (tol^2) * gg0
            info = 0;
            break;
        end
        sg = s'*g;
        shs = s'*hs;
        sgk = sg + shs;
        if sgk / sqrt(gg * delsq) <= tol - 1
            info = 0;
            break;
        end

        % Begin the 2D minimization by calculating D and HD and some
        % scalar products.
        t = sqrt(delsq * gg - sgk * sgk);
        d = (delsq / t) * (g + hs) - (sgk / t) * s;
        %----------------------------------------------------------------%
        %-----%hd = hessmul(hq, pq, xpt, d) %----------------------------%
        hd = xpt * (pq .* (d' * xpt)') + hq * d;
        %----------------------------------------------------------------%
        dg = d'*g;
        dhd = hd'*d;
        dhs = hd'*s;

        % Seek the value of the angle that minimizes Q.
        cf = 0.5 * (shs - dhd);
        qbeg = sg + cf;
        qsav = qbeg;
        qmin = qbeg;
        isav = 0;
        iu = 49;
        unitang = (2 * pi) / (iu + 1);

        for i = 1:iu
            angle = i * unitang;
            cth = cos(angle);
            sth = sin(angle);
            qnew = (sg + cf * cth) * cth + (dg + dhs * cth) * sth;
            if qnew < qmin
                qmin = qnew;
                isav = i;
                quada = qsav;
            elseif i == isav + 1
                quadb = qnew;
            end
            qsav = qnew;
        end

        if isav == 0
            quada = qnew;
        end
        if isav == iu
            quadb = qbeg;
        end
        if abs(quada - quadb) > 0
            quada = quada - qmin;
            quadb = quadb - qmin;
            angle = 0.5 * (quada - quadb) / (quada + quadb);
        else
            angle = 0;
        end
        angle = unitang * (isav + angle);

        % Calculate the new S.
        cth = cos(angle);
        sth = sin(angle);
        reduc = qbeg - (sg + cf * cth) * cth - (dg + dhs * cth) * sth;
        sold = s;
        s = cth * s + sth * d;

        % Exit in case of Inf/NaN in S.
        if isinf(sum(abs(s)))
            s = sold;
            info = -1;
            break;
        end

        % Calculate HS. Then test for convergence.
        hs = cth * hs + sth * hd;
        gg = norm(g + hs)^2;
        qred = qred + reduc;
        if reduc / qred <= tol
            info = 1;
            break;
        end
    end

end
