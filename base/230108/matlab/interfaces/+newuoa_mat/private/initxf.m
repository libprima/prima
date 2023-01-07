function [ij, kopt, nf, fhist, fval, xbase, xhist, xpt, info] = initxf(calfun, iprint, ftarget, rhobeg, x, maxfhist, maxxhist, debugflag, npt)
    n = length(x);
    % INITXF performs the initialization regarding the interpolation
    % points and corresponding function values.

    % Solver-specific module

    ij = NaN(2, npt);
    fval = NaN(npt, 1);
    fhist = NaN(1, maxfhist);
    xbase = NaN(n, 1);
    xhist = NaN(n, maxxhist);
    xpt = NaN(n, npt);
    % Remark on IJ:
    % When K > 2*N + 1, all the entries of XPT(:, K) will be zero except
    % that the IJ(1, K) and IJ(2, K) entries will be set to RHOBEG or
    % -RHOBEG. Consequently, the Hessian of the quadratic model will get a
    % possibly nonzero (IJ(1, K), IJ(2, K)) entry.
    % Indeed, IJ(:, 1 : 2*N + 1) is never used.

    % Local variables
    solver = 'NEWUOA';
    funname = 'INITXF';

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(X) or SIZE(FVAL) is invalid.', funname);
        end
        if size(xhist, 1) ~= n && maxxhist > 0
            error('Error: %s: XHIST is nonempty but SIZE(XHIST, 1) /= SIZE(X).', funname);
        end
        if maxfhist * maxxhist > 0 && maxfhist ~= maxxhist
            error('Error: %s: FHIST and XHIST are both nonempty but SIZE(FHIST) /= SIZE(XHIST, 2).', funname);
        end
        verisize(ij, 2, npt);
        verisize(xbase, n);
        verisize(xpt, n, npt);
    end

    % At return,
    % INFO = 0: initialization finishes normally
    % INFO = infos('ftarget_achieved'): return because f <= ftarget
    % INFO = infos('nan_x'): return because x contains NaN
    % INFO = infos('nan_inf_f'): return because f is either NaN or +infinity
    info = 0;

    % We set ij = 1 in case the initialization aborts due to abnormality. If
    % we do not do this, ij will be undefined if the initialization aborts.
    ij = ones(2, npt);

    % Set XBASE to X.
    xbase = x;

    % Initialize XPT to 0.
    xpt = zeros(n, npt);

    % Begin the initialization procedure. The coordinates of the
    % displacement of the next initial interpolation point from XBASE are
    % set in XPT(:, .).

    % EVALUATED is a boolean array indicating whether the function value of
    % the i-th interpolation point has been evaluated. We need it for a
    % portable counting of the number of function evaluations, especially
    % if the loop is conducted asynchronously. However, the loop here is
    % not fully parallelizable if NPT>2N+1, because the definition
    % XPT(;, 2N+2:end) depends on FVAL(1:2N+1).
    evaluated = false(length(fval), 1);

    % NPT_REVISED equals NPT, unless it turns out necessary to return due to
    % abnormality (NaN or Inf occurs, or F < FTARGET).
    npt_revised = npt;

    % Set XPT, FVAL, KOPT, and XOPT.

    % Set XPT(:, 2 : N + 1).
    for k = 2:min(npt, n + 1)
        xpt(k - 1, k) = rhobeg;
    end
    % Set XPT(:, N+2 : NPT)
    for k = (n + 2):min(npt, 2 * n + 1)
        xpt(k - n - 1, k) = -rhobeg;
    end

    % Set FVAL(1 : 2*N + 1) by evaluating F. Totally parallelizable except for
    % FMSSG, which outputs messages to the console or files.
    for k = 1:min(npt, 2 * n + 1)
        xtemp = xpt(:, k) + xbase;
        if any(isnan(xtemp))
            info = infos('nan_x');
            npt_revised = 0;
            break;
        end
        f = calfun(xtemp);
        fmssg(iprint, k, f, xtemp, solver);
        evaluated(k) = true;
        fval(k) = f;

        if maxfhist >= 1
            khist = mod(k - 1, maxfhist) + 1;
            fhist(khist) = f;
        end
        if maxxhist >= 1
            khist = mod(k - 1, maxxhist) + 1;
            xhist(:, khist) = xtemp;
        end

        % Check whether to exit.
        if f <= ftarget
            info = infos('ftarget_achieved');
            npt_revised = 0;
            break;
        end
        if (f == Inf) || isnan(f)
            info = infos('nan_inf_f');
            npt_revised = 0;
            break;
        end
    end

    % Set XPT(:, 2*N + 2 : NPT). It depends on FVAL(2 : 2*N + 1).
    for k = (2 * n + 2):npt_revised
        % Decide IJ(:, K).  In general, when NPT = (N+1)*(N+2)/2, we can
        % set IJ(1 : NPT - (2*N + 1)) to ANY permutation of
        % {(I, J) : 1 <= I > J <= N};
        % when NPT < (N+1)*(N+2)/2, IJ(1 : NPT - (2*N + 1)) is the
        % first NPT - (2*N - 1) elements of such a permutation. Powell took
        % the following permutation:
        itemp = (k - n - 2) / n;
        j = k - (itemp + 1) * n - 1;
        i = j + itemp;
        if i > n
            itemp = j;
            j = i - n;
            i = itemp;
        end
        ij(1, k) = i;
        ij(2, k) = j;

        % The following lines set XPT(;, K) to
        % XPT(:, I + 1) or XPT(:, I + N + 1)
        % +
        % XPT(:, J + 1) or XPT(:, J + N + 1),
        % depending on the values of FVAL(I + 1), FVAL(I + N + 1),
        % FVAL(J + 1), and FVAL(J + N + 1).
        %
        % This is the only places where the definition
        % of XPT(:, 2*N + 2 : NPT) depends on F(2 : 2*N + 1).
        % If we set XPT(:, K) to XPT(:, I + 1) + XPT(:, J + 1)
        % regardless of FVAL, then the evaluations of FVAL(1 : NPT)
        % can be merged, and they are totally parallelizable; this can be
        % benificial if the function evaluations are expensive, which is
        % likely the case.
        if fval(i + 1) <= fval(i + n + 1)
            xpt(i, k) = xpt(i, i + 1);
        else
            xpt(i, k) = xpt(i, i + n + 1);
        end
        if fval(j + 1) <= fval(j + n + 1)
            xpt(j, k) = xpt(j, j + 1);
        else
            xpt(j, k) = xpt(j, j + n + 1);
        end
    end

    % Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable except
    % FMSSG, which outputs messages to the console or files.
    for k = (2 * n + 2):npt_revised
        xtemp = xpt(:, k) + xbase;
        if any(isnan(xtemp))
            info = infos('nan_x');
            break;
        end
        f = calfun(xtemp);
        fmssg(iprint, k, f, xtemp, solver);
        evaluated(k) = true;
        fval(k) = f;

        if maxfhist >= 1
            khist = mod(k - 1, maxfhist) + 1;
            fhist(khist) = f;
        end
        if maxxhist >= 1
            khist = mod(k - 1, maxxhist) + 1;
            xhist(:, khist) = xtemp;
        end

        % Check whether to exit.
        if f <= ftarget
            info = infos('ftarget_achieved');
            break;
        end
        if (f == Inf) || isnan(f)
            info = infos('nan_inf_f');
            break;
        end
    end

    % Set NF, KOPT
    nf = sum(evaluated);
    minimum = min(fval(evaluated));
    kopt = find(fval<= minimum, 1);

end
