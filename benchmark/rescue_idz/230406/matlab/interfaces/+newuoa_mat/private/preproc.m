function [iprint, maxfun, maxhist, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend] = preproc(n, iprint, maxfun, maxhist, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend)

    % Local variable
    solver = 'NEWUOA';

    if iprint ~= 0 && abs(iprint) ~= 1 && abs(iprint) ~= 2 && abs(iprint) ~= 3
        iprint = consts('iprint_dft');
        fprintf('%s: invalid IPRINT; it should be 0, 1, -1, 2, -2, 3, or -3; it is set to %d.\n', solver, iprint);
    end

    if maxfun < n + 3
        maxfun = n + 3;
        fprintf('%s: invalid MAXFUN; it should an integer at least N + 3 ; it is set to %d.\n', solver, maxfun);
    end

    if maxhist < 0
        maxhist = maxfun;
        fprintf('%s: invalid MAXHIST; it should be a nonnegative integer; it is set to %d.\n', solver, maxhist);
    end
    % MAXHIST > MAXFUN is never needed.
    maxhist = min(maxhist, maxfun);

    if npt < n + 2 || npt > min(maxfun - 1, ((n + 2) * (n + 1)) / 2)
        npt = min(maxfun - 1, 2 * n + 1);
        fprintf('%s: invalid NPT; it should an integer in the interval [N+2, (N+1)(N+2)/2], and it should be less than MAXFUN; it is set to %d.\n', solver, npt);
    end

    if isnan(ftarget)
        ftarget = consts('ftarget_dft');
        fprintf('%s: invalid FTARGET; it should a real number; it is set to %15.6f.\n', solver, ftarget);
    end

    % When the difference between ETA1 and ETA2 is tiny, we force them to equal.
    % See the explanation around RHOBEG and RHOEND for the reason.
    if abs(eta1 - eta2) < 1.0e2 * eps * max(abs(eta1), 1)
        eta2 = eta1;
    end

    if isnan(eta1)
        % In this case, we take the value hard coded in Powell's original code
        % without any warning. It is useful when interfacing with MATLAB/Python.
        eta1 = 0.1;
    elseif eta1 < 0 || eta1 >= 1
        % Take ETA1 into account if it has a valid value.
        if eta2 > 0 && eta2 <= 1
            eta1 = max(eps, eta2 / 7.0);
        else
            eta1 = 0.1;
        end
        fprintf('%s: invalid ETA1; it should be in the interval [0, 1) and not more than ETA2; it is set to %15.6f.\n', solver, eta1);
    end

    if isnan(eta2)
        % In this case, we take the value hard coded in Powell's original code
        % without any warning. It is useful when interfacing with MATLAB/Python.
        eta2 = 0.7;
    elseif eta2 < eta1 || eta2 > 1
        eta2 = (eta1 + 2) / 3.0;
        fprintf('%s: invalid ETA2; it should be in the interval [0, 1] and not less than ETA1; it is set to %15.6f.\n', solver, eta2);
    end

    if isnan(gamma1)
        % In this case, we take the value hard coded in Powell's original code
        % without any warning. It is useful when interfacing with MATLAB/Python.
        gamma1 = 0.5;
    elseif gamma1 <= 0 || gamma1 >= 1
        gamma1 = 0.5;
        fprintf('%s: invalid GAMMA1; it should in the interval (0, 1); it is set to %15.6f.\n', solver, gamma1);
    end

    if isnan(gamma2)
        % In this case, we take the value hard coded in Powell's original code
        % without any warning. It is useful when interfacing with MATLAB/Python.
        gamma2 = 2;
    elseif gamma2 < 1 || isinf(gamma2)
        gamma2 = 2;
        fprintf('%s: invalid GAMMA2; it should a real number not less than 1; it is set to %15.6f.\n', solver, gamma2);
    end

    if abs(rhobeg - rhoend) < 1.0e2 * eps * max(abs(rhobeg), 1)
        % When the data is passed from the interfaces (e.g., MEX) to the Fortran
        % code, RHOBEG, and RHOEND may change a bit. It was observed in a MATLAB
        % test that MEX passed 1 to Fortran as 0.99999999999999978. Therefore,
        % if we set RHOEND = RHOBEG in the interfaces, then it may happen that
        % RHOEND > RHOBEG, which is considered as an invalid input. To avoid this
        % situation, we force RHOBEG and RHOEND to equal when the difference is tiny.
        rhoend = rhobeg;
    end

    if rhobeg <= 0 || isnan(rhobeg) || isinf(rhobeg)
        % Take RHOEND into account if it has a valid value.
        if ~isinf(rhoend) && rhoend > 0
            rhobeg = max(10 * rhoend, consts('rhobeg_dft'));
        else
            rhobeg = consts('rhobeg_dft');
        end
        fprintf('%s: invalid RHOBEG; it should be a positive number; it is set to %15.6f.\n', solver, rhobeg);
    end

    if rhoend <= 0 || rhobeg < rhoend || isnan(rhoend) || isinf(rhoend)
        rhoend = max(eps, min(0.1 * rhobeg, consts('rhoend_dft')));
        fprintf('%s: invalid RHOEND; it should be a positive number and RHOEND <= RHOBEG; it is set to %15.6f.\n', solver, rhoend);
    end

end
