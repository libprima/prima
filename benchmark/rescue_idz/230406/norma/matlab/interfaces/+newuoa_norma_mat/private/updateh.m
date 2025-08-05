function [idz, bmat, zmat] = updateh(knew, beta, vlag_in, idz, bmat, zmat, debugflag)
    n = size(bmat, 1);
    npt = size(bmat, 2) - n;
    % UPDATE updates arrays BMAT and ZMAT together with IDZ, in order to replace the interpolation point
    % XPT(:, KNEW) by XNEW. On entry, VLAG_IN contains the components of the vector THETA*WCHECK + e_b
    % of the updating formula (6.11) in the NEWUOA paper, and BETA holds the value of the parameter that
    % has this name. VLAG_IN and BETA contains information about XNEW, because they are calculated
    % according to D = XNEW - XOPT.
    %
    % See Section 4 of the NEWUOA paper.

    % Local variables
    w = NaN(length(vlag_in), 1);
    funname = 'UPDATEH';

    % Get and verify the sizes.

    if debugflag
        if n == 0 || npt < n + 2
            error('Error: %s: SIZE(BMAT) is invalid.', funname);
        end
        verisize(zmat, npt, npt - n - 1);
        verisize(vlag_in, npt + n);
    end

    vlag = vlag_in;    % VLAG_IN is IN10T(IN) and cannot be revised.

    % Apply rotations to put zeros in the KNEW-th row of ZMAT. A 2x2 rotation will be multiplied to ZMAT
    % from the right so that ZMAT(KNEW, [JL, J]) becomes [SQRT(ZMAT(KNEW, JL)^2 + ZMAT(KNEW, J)^2), 0].
    % As in MATLAB, PLANEROT(X) returns a 2x2 Givens matrix G for X in R^2 so that Y = G*X has Y(2) = 0.

    % In the loop, if 2 <= J < IDZ, then JL = 1; if IDZ < J <= NPT - N - 1, then JL = IDZ.
    jl = 1;
    for j = 2:(npt - n - 1)
        if j == idz
            jl = idz;
            break;
        end
        if abs(zmat(knew, j)) > 0
            grot = planerot(zmat(knew, [jl; j])');            % MATLAB code: GROT = PLANEROT(ZMAT(KNEW, [JL, J])')
            zmat(:, [jl; j]) = zmat(:, [jl; j]) * grot';
            zmat(knew, j) = 0;
        end
    end
    % The value of JL after the loop is important below. Its value is determined by the current (i.e.,
    % unupdated) value of IDZ. IDZ is an integer in {1, ..., NPT-N} such that s_j = -1 for j < IDZ while
    % s_j = 1 for j >= IDZ in the factorization of OMEGA. See (3.17) and (4.16) of the NEWUOA paper.
    % The value of JL has two possibilities:
    % 1. JL = 1 iff IDZ = 1 or IDZ = NPT - N.
    % 1.1. IDZ = 1 means that
    % OMEGA = sum_{J=1}^{NPT-N-1} ZMAT(:, J)*ZMAT(:, J)' ;
    % 1.2. IDZ = NPT - N means that
    % OMEGA = - sum_{J=1}^{NPT-N-1} ZMAT(:, J)*ZMAT(:, J)' ;
    % 2. JL = IDZ > 1 iff 2 <= IDZ <= NPT - N - 1.

    % Put the first NPT components of the KNEW-th column of HLAG into W, and calculate the parameters of
    % the updating formula.
    tempa = zmat(knew, 1);
    if idz >= 2
        tempa = -tempa;
    end

    w(1:npt) = tempa * zmat(:, 1);
    if jl > 1
        tempb = zmat(knew, jl);
        w(1:npt) = w(1:npt) + tempb * zmat(:, jl);
    end

    alpha = w(knew);
    tau = vlag(knew);
    tausq = tau * tau;
    denom = alpha * beta + tausq;
    % After the following line, VLAG = Hw - e_t in the NEWUOA paper.
    vlag(knew) = vlag(knew) - 1;
    sqrtdn = sqrt(abs(denom));

    % Complete the updating of ZMAT when there is only one nonzero element in the KNEW-th row of the new
    % matrix ZMAT, but, if IFLAG is set to one, then the first column of ZMAT will be exchanged with
    % another one later.
    reduce_idz = false;
    if jl == 1
        % There is only one nonzero in ZMAT(KNEW, :) after the rotation. This is the normal case,
        % because IDZ is 1 in precise arithmetic.
        %---------------------------------------------------------------------------------------------%
        % Up to now, TEMPA = ZMAT(KNEW, 1) if IDZ = 1 and TEMPA = -ZMAT(KNEW, 1) if IDZ >= 2. However,
        % according to (4.18) of the NEWUOA paper, TEMPB should always be ZMAT(KNEW, 1)/sqrtdn
        % regardless of IDZ. Therefore, the following definition of TEMPB is inconsistent with (4.18).
        % This is probably a BUG. See also Lemma 4 and (5.13) of Powell's paper "On updating the inverse
        % of a KKT matrix". However, the inconsistency is hardly observable in practice, because JL = 1
        % implies IDZ = 1 in precise arithmetic.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
        %tempb = tempa/sqrtdn
        %tempa = tau/sqrtdn
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
        % Here is the corrected version (only TEMPB is changed).
        tempa = tau / sqrtdn;
        tempb = zmat(knew, 1) / sqrtdn;
        %---------------------------------------------------------------------------------------------%

        zmat(:, 1) = tempa * zmat(:, 1) - tempb * vlag(1:npt);

        %---------------------------------------------------------------------------------------------%
        % The following six lines by Powell are obviously problematic --- TEMP is always nonnegative.
        % According to (4.18) of the NEWUOA paper, the "TEMP < 0" and "TEMP >= 0" below should be
        % both revised to "DENOM < 0". See also the corresponding part of the LINCOA code. Note that
        % the NEAUOA paper uses SIGMA to denote DENOM. Check also Lemma 4 and (5.13) of Powell's paper
        % "On updating the inverse of a KKT matrix". It seems that the BOBYQA code does not have this
        % part --- it does not have IDZ at all (why?). Anyway, these lines are not invoked very often in
        % practice, because IDZ should always be 1 in precise arithmetic.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
        %if (idz == 1 .and. sqrtdn < 0) then
        %    idz = 2
        %end if
        %if (idz >= 2 .and. sqrtdn >= 0) then
        %    reduce_idz = .true.
        %end if
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
        % This is the corrected version. It copies precisely the
        % corresponding part of the LINCOA code.
        if denom < 0
            if idz == 1
                % This is the first place (out of two) where IDZ is
                % increased. Note that IDZ = 2 <= NPT-N after the update.
                idz = 2;
            else
                % This is the first place (out of two) where IDZ is
                % decreased (by 1). Since IDZ >= 2 in this case, we have
                % IDZ >= 1 after the update.
                reduce_idz = true;
            end
        end
        %---------------------------------------------------------------------------------------------%
    else
        % Complete the updating of ZMAT in the alternative case. ZMAT(KNEW, :) has two nonzeros after
        % the rotations.
        ja = 1;
        if beta >= 0
            ja = jl;
        end
        jb = jl + 1 - ja;
        temp = zmat(knew, jb) / denom;
        tempa = temp * beta;
        tempb = temp * tau;
        temp = zmat(knew, ja);
        scala = 1 / sqrt(abs(beta) * temp * temp + tausq);
        scalb = scala * sqrtdn;
        zmat(:, ja) = scala * (tau * zmat(:, ja) - temp * vlag(1:npt));
        zmat(:, jb) = scalb * (zmat(:, jb) - tempa * w(1:npt) - tempb * vlag(1:npt));
        % If and only if DENOM <= 0, IDZ will be revised according to the sign of BETA.
        % See (4.19)--(4.20) of the NEWUOA paper.
        if denom <= 0
            if beta < 0
                % This is the second place (out of two) where IDZ is increased. Since
                % JL = IDZ <= NPT-N-1 in this case, we have IDZ <= NPT-N after the update.
                idz = idz + 1;
            end
            if beta >= 0
                % This is the second place (out of two) where IDZ is decreased (by 1). Since IDZ >= 2
                % in this case, we have IDZ >= 1 after the update.
                reduce_idz = true;
            end
        end
    end

    % IDZ is reduced in the following case, and usually the first column of ZMAT is exchanged with a
    % later one.
    if reduce_idz
        idz = idz - 1;
        if idz > 1
            ztemp = zmat(:, 1);
            zmat(:, 1) = zmat(:, idz);
            zmat(:, idz) = ztemp;
        end
    end

    % Finally, update the matrix BMAT.
    w(npt + 1:npt + n) = bmat(:, knew);
    v1 = (alpha * vlag(npt + 1:npt + n) - tau * w(npt + 1:npt + n)) / denom;
    v2 = (-beta * w(npt + 1:npt + n) - tau * vlag(npt + 1:npt + n)) / denom;

    bmat = bmat + (1 * (v1 * vlag') + 1 * (v2 * w'));
    % In floating-point arithmetic, the update above does not guarantee BMAT(:, NPT+1 : NPT+N) to be
    % symmetric. Symmetrization needed.
    bmat(:, npt + 1:npt + n) = (bmat(:, npt + 1:npt + n)' + bmat(:, npt + 1:npt + n)) / 2;

end
