function [x, nf, f, fhist, xhist, info] = newuob(calfun, iprint, maxfun, npt, eta1, eta2, ftarget, gamma1, gamma2, rhobeg, rhoend, x, maxfhist, maxxhist, debugflag)
    n = length(x);
    xhist = NaN(n, maxxhist);
    % NEWUOB performs the actual calculations of NEWUOA. The arguments IPRINT, MAXFUN, MAXHIST, NPT,
    % ETA1, ETA2, FTARGET, GAMMA1, GAMMA2, RHOBEG, RHOEND, X, NF, F, FHIST, XHIST, and INFO are
    % identical to the corresponding arguments in subroutine NEWUOA.

    % XBASE holds a shift of origin that should reduce the contributions from rounding errors to values
    % of the model and Lagrange functions.
    % XOPT is the displacement from XBASE of the best vector of variables so far (i.e., the one provides
    % the least calculated F so far). FOPT = F(XOPT + XBASE).
    % D is reserved for trial steps from XOPT.
    % XNEW = XOPT + D, corresponding to the vector of variables for the next calculation of F.
    % [XPT, FVAL, KOPT] describes the interpolation set:
    % XPT contains the interpolation points relative to XBASE, each COLUMN for a point; FVAL holds the
    % However, there is a delay between the update of XOPT and KOPT. So they are not always consistent
    % in the mid of an iteration. See the comment on the update of XOPT for details.
    % values of F at the interpolation points; KOPT is the index of XOPT in XPT (XPT(:,KOPT) = XOPT).
    % [GQ, HQ, PQ] describes the quadratic model: GQ will hold the gradient of the quadratic model at
    % XBASE; HQ will hold the explicit second order derivatives of the quadratic model; PQ will contain
    % the parameters of the implicit second order derivatives of the quadratic model.
    % [BMAT, ZMAT, IDZ] describes the matrix H in the NEWUOA paper (eq. 3.12), which is the inverse of
    % the coefficient matrix of the KKT system for the least-Frobenius norm interpolation problem:
    % BMAT will hold the last N ROWs of H; ZMAT will hold the factorization of the leading NPT*NPT
    % submatrix of H, this factorization being ZMAT*Diag(DZ)*ZMAT^T with DZ(1:IDZ-1)=-1, DZ(IDZ:NPT)=1.
    % VLAG will contain the values of the Lagrange functions at a new point X. They are part of a
    % product that requires VLAG to be of length NPT + N. Both VLAG and BETA are critical for the
    % updating procedure of H, which is detailed formula (4.10)--(4.12) of the NEWUOA paper.
    %
    % See Section 2 of the NEWUOA paper for more information about these variables.

    % Local variables
    d = NaN(length(x), 1);
    dnormsav = NaN(3, 1);
    solver = 'NEWUOA';
    funname = 'NEWUOB';

    % Get and verify the sizes.

    if debugflag
        if n == 0
            error('Error: %s: SIZE(X) is invalid.', funname);
        end
        if size(xhist, 1) ~= n && maxxhist > 0
            error('Error: %s: XHIST is nonempty but SIZE(XHIST, 1) /= SIZE(X).', funname);
        end
        if maxfhist * maxxhist > 0 && maxfhist ~= maxxhist
            error('Error: %s: FHIST and XHIST are both nonempty but SIZE(FHIST) /= SIZE(XHIST, 2).', funname);
        end
    end

    % Initialize FVAL, XBASE, and XPT.
    [ij, kopt, nf, fhist, fval, xbase, xhist, xpt, subinfo] = initxf(calfun, iprint, ftarget, rhobeg, x, maxfhist, maxxhist, debugflag, npt);
    xopt = xpt(:, kopt);
    fopt = fval(kopt);
    x = xbase + xopt;    % Set X.
    f = fopt;    % Set F.

    % Check whether to return after initialization.
    if subinfo == infos('ftarget_achieved') || subinfo == infos('nan_x') || subinfo == infos('nan_inf_f')
        % In these cases, pack the data and return immediately.
        info = subinfo;
        % Rearrange FHIST and XHIST so that they are in the chronological order.
        if maxfhist >= 1 && maxfhist < nf
            khist = mod(nf - 1, maxfhist) + 1;
            fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)];
        end
        if maxxhist >= 1 && maxxhist < nf
            khist = mod(nf - 1, maxxhist) + 1;
            xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)];
            % The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
            % order of Fortran arrays. In MATLAB, it is NOT necessary to call `reshape`.
        end
        retmssg(info, iprint, nf, f, x, solver);
        return;
    end

    % Initialize GQ, HQ, and PQ.
    [gq, hq, pq] = initq(ij, fval, xpt, debugflag);

    % Initialize BMAT and ZMAT, and IDZ.
    [idz, bmat, zmat] = inith(ij, xpt, debugflag);

    % After initializing GQ, HQ, PQ, BMAT, ZMAT, one can also choose to return if subinfo = infos('nan_model')
    % (NaN occurs in the model). We do not do it here. If such a model is harmful, then it will probably
    % lead to other returns (NaN in X, NaN in F, trust region subproblem fails, ...); otherwise, the
    % code will continue to run and possibly get rid of the NaN in the model.

    % Set some more initial values and parameters.
    rho = rhobeg;
    delta = rho;
    moderrsav = Inf(length(dnormsav), 1);
    dnormsav = Inf(3, 1);
    itest = 0;
    trtol = 1.0e-2;    % Tolerance used in trsapp.
    % We must initialize RATIO. Otherwise, when SHORTD = TRUE, compilers will raise a run-time error
    % that RATIO is undefined. Powell's code indeed sets RATIO = -1 when SHORTD is TRUE, and use this
    % artificial value when setting IMPROVE_GEO and REDUCE_RHO; however, we choose not to use the
    % artificial RATIO but use SHORTD to be more explicit. See IMPROVE_GEO and REDUCE_RHO for details.
    ratio = -1;

    % Normally, each trust-region iteration takes one function evaluation. The following setting
    % essentially imposes no constraint on the maximal number of trust region iterations.
    maxtr = 10 * maxfun;
    % MAXTR is unlikely to be reached, but we define the following default value for INFO for safety.
    info = infos('maxtr_reached');

    % Begin the iterative procedure.
    % After solving a trust-region subproblem, NEWUOA uses 3 boolean variables to control the work flow.
    % SHORTD - Is the trust region trial step too short to invoke a function evaluation?
    % IMPROVE_GEO - Will we improve the model after the trust region iteration?
    % REDUCE_RHO - Will we reduce rho after the trust region iteration?
    % REDUCE_RHO = REDUCE_RHO_1 .OR. REDUCE_RHO_2 (see boxes 14 and 10 of Fig. 1 in the NEWUOA paper).
    % NEWUOA never sets IMPROVE_GEO and REDUCE_RHO to TRUE simultaneously.
    for tr = 1:maxtr
        % Solve the trust region subproblem.
        [crvmin, d] = trsapp(delta, gq, hq, pq, trtol, xopt, xpt, debugflag);

        % Calculate the length of the trial step D.
        dnorm = min(delta, norm(d));

        % SHORTD corresponds to Box 3 of the NEWUOA paper.
        shortd = (dnorm < 0.5 * rho);
        % REDUCE_RHO_1 corresponds to Box 14 of the NEWUOA paper.
        reduce_rho_1 = shortd && (max(abs(moderrsav)) <= 0.125 * crvmin * rho * rho) && (max(dnormsav) <= rho);
        if shortd && (~reduce_rho_1)
            % Reduce DELTA. After this, DELTA < DNORM may hold.
            delta = 0.1 * delta;
            if delta <= 1.5 * rho
                delta = rho;                % Set DELTA to RHO when it is close.
            end
        end

        if ~shortd            % D is long enough.
            % Shift XBASE if XOPT may be too far from XBASE.
            %if (inprod(d, d) <= 1.0e-3_RP*inprod(xopt, xopt)) then  % Powell
            if dnorm * dnorm <= 1.0e-3 * norm(xopt)^2
                [bmat, gq, hq, xbase, xopt, xpt] = shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt, debugflag);
            end

            % Calculate VLAG and BETA for D. It makes uses of XOPT, so this is done before updating XOPT.
            [beta, vlag] = vlagbeta(idz, kopt, bmat, d, xpt, zmat, debugflag);

            % Use the current quadratic model to predict the change in F due to the step D.
            qred = calquad(d, gq, hq, pq, xopt, xpt, debugflag);

            % Calculate the next value of the objective function.
            xnew = xopt + d;
            x = xbase + xnew;
            if any(isnan(x))
                f = sum(x);                % Set F to NaN. It is necessary.
                info = infos('nan_x');
                break;
            end
            f = calfun(x);
            nf = nf + 1;
            fmssg(iprint, nf, f, x, solver);
            if maxfhist >= 1
                khist = mod(nf - 1, maxfhist) + 1;
                fhist(khist) = f;
            end
            if maxxhist >= 1
                khist = mod(nf - 1, maxxhist) + 1;
                xhist(:, khist) = x;
            end

            % DNORMSAVE contains the DNORM of the latest 3 function evaluations with the current RHO.
            dnormsav = [dnormsav(2:length(dnormsav)); dnorm];

            % MODERR is the error of the current model in predicting the change in F due to D.
            moderr = f - fopt - qred;
            % MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
            moderrsav = [moderrsav(2:length(moderrsav)); moderr];

            % Calculate the reduction ratio and update DELTA accordingly.
            if isnan(qred) || qred >= 0
                info = infos('trsubp_failed');
                break;
            end
            ratio = (f - fopt) / qred;
            % Update DELTA. After this, DELTA < DNORM may hold.
            delta = trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio);
            if delta <= 1.5 * rho
                delta = rho;
            end

            % Update XOPT and FOPT. Before KOPT is updated, XOPT may differ from XPT(:, KOPT), and FOPT
            % may differ from FVAL(KOPT). Note that the code may exit before KOPT is updated. See below.
            % The updated XOPT is needed by SETDROP_TR.
            if f < fopt
                xopt = xnew;
                fopt = f;
            end

            % Check whether to exit
            if isnan(f) || (f == Inf)
                info = infos('nan_inf_f');
                break;
            end
            if f <= ftarget
                info = infos('ftarget_achieved');
                break;
            end
            if nf >= maxfun
                info = infos('maxfun_reached');
                break;
            end

            % Set KNEW_TR to the index of the interpolation point that will be replaced by XNEW. KNEW_TR
            % will ensure that the geometry of XPT is "good enough" after the replacement. Note that the
            % information of XNEW is included in VLAG and BETA, which are calculated according to
            % D = XNEW - XOPT. KNEW_TR = 0 means it is impossible to obtain a good interpolation set
            % by replacing any current interpolation point with XNEW.
            knew_tr = setdrop_tr(idz, kopt, beta, delta, ratio, rho, vlag(1:npt), xopt, xpt, zmat, debugflag);

            if knew_tr > 0
                % If KNEW_TR > 0, then update BMAT, ZMAT and IDZ, so that the KNEW_TR-th interpolation
                % point is replaced by XNEW. If KNEW_TR = 0, then probably the geometry of XPT needs
                % improvement, which will be handled below.
                [idz, bmat, zmat] = updateh(knew_tr, beta, vlag, idz, bmat, zmat, debugflag);

                % Update the quadratic model using the updated BMAT, ZMAT, IDZ.
                [gq, hq, pq] = updateq(idz, knew_tr, bmat(:, knew_tr), moderr, zmat, xpt(:, knew_tr), gq, hq, pq, debugflag);

                % Include the new interpolation point.
                xpt(:, knew_tr) = xnew;                % Should be done after UPDATEQ.
                fval(knew_tr) = f;
                if fval(knew_tr) < fval(kopt)
                    kopt = knew_tr;
                end
                % KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW_TR) = FVAL(KOPT) and
                % KNEW_TR < KOPT, then MINLOC(FVAL) = KNEW_TR /= KOPT. Do not change KOPT in this case.
            end

            % Test whether to replace the new quadratic model Q by the least-Frobenius norm interpolant
            % Q_alt. Perform the replacement if certain criteria are satisfied.
            % N.B.:
            % 1. This part is OPTIONAL, but it is crucial for the performance on some problems. See
            % Section 8 of the NEWUOA paper.
            % 2. TRYQALT is called only after a trust-region step but not after a geometry step, maybe
            % because the model is expected to be good after a geometry step.
            % 3. If KNEW_TR = 0 after a trust-region step, TRYQALT is not invoked. In this case, the
            % interpolation set is unchanged, so it seems reasonable to keep the model unchanged.
            % 4. In theory, FVAL - FOPT in the call of TRYQALT can be replaced by FVAL + C with any
            % constant C. This constant will not affect the result in precise arithmetic. Powell chose
            % C = - FVAL(KOPT_OLD), where KOPT_OLD is the KOPT before the update above (Powell updated
            % KOPT after TRYQALT). Here we use C = -FOPT, as it worked slightly better on CUTEst,
            % although there is no difference theoretically. Note that FVAL(KOPT_OLD) may not equal
            % FOPT_OLD --- it may happen that KNEW_TR = KOPT_OLD so that FVAL(KOPT_OLD) has been revised
            % after the last function evaluation.
            % 5. Question: Since TRYQALT is invoked only when DELTA equals the current RHO, why not
            % reset ITEST to 0 when RHO is reduced?
            if knew_tr > 0 && delta <= rho                % DELTA = RHO.
                [itest, gq, hq, pq] = tryqalt(idz, fval - fopt, ratio, bmat(:, 1:npt), zmat, itest, gq, hq, pq, debugflag);
            end
        end        % End of if (.not. shortd)

        % Before next trust region iteration, we may improve the geometry of XPT or reduce rho
        % according to IMPROVE_GEO and REDUCE_RHO. Now we decide these two indicators.

        % First define IMPROVE_GEO, which corresponds to Box 8 of the NEWUOA paper.
        % The geometry of XPT likely needs improvement if the trust-region step bad --- either too short
        % (SHORTD = TRUE) or the reduction ratio is small (RATIO < 0.1_RP). However, if REDUCE_RHO_1 is
        % TRUE, meaning that the step is short and the latest model errors have been small, then we do
        % not need to improve the geometry; instead, RHO will be reduced.
        % To improve the geometry of XPT, we will check whether the interpolation points are all close
        % enough to the best point so far, i.e., all the points are within a ball centered at XOPT with
        % a radius of 2*DELTA. If not, the farthest point will be replaced with a point selected by
        % GEOSTEP, aiming to ameliorate the geometry of the interpolation set; if yes, then RHO will be
        % reduced if MAX(DELTA, DNORM) <= RHO (if MAX(DELTA, DNORM) > RHO, then, as Powell mentioned
        % under (2.3) of the NEWUOA paper, "RHO has not restricted the most recent choice of D", so it
        % is not reasonable to reduce RHO).
        % N.B.:
        % 1. RATIO must be set even if SHORTD = TRUE. Otherwise, compilers will raise a run-time error.
        % 2. If SHORTD = FALSE and KNEW_TR = 0, then IMPROVE_GEO = TRUE. Therefore, IMPROVE_GEO = TRUE
        % if it is impossible to obtain a good XPT by replacing a current point with the one suggested
        % by the trust region step.
        % 3. If REDUCE_RHO = FALSE and SHORTD = TRUE, then the trust-region step is not tried at all,
        % i.e., no function evaluation is invoked at XOPT + D (when REDUCE_RHO = TRUE, the step is not
        % tried either, but the same step will be generated again at the next trust-region iteration
        % after RHO is reduced and DELTA is updated; see the end of Section 2 of the NEWUOA paper).
        % 4. If SHORTD = FALSE and KNEW_TR = 0, then the trust-region step invokes a function evaluation
        % at XOPT + D, but [XOPT + D, F(XOPT +D)] is not included into [XPT, FVAL]. In other words, this
        % function value is discarded. THEORETICALLY, KNEW_TR = 0 only if RATIO <= 0, so that a function
        % value that renders a reduction is never discarded; however, KNEW_TR may turn out 0 due to NaN
        % even if RATIO > 0. See SETDROP_TR for details.
        % 5. If SHORTD = FALSE and KNEW_TR > 0 and RATIO < 0.1_RP, then [XPT, FVAL] is updated so that
        % [XPT(KNEW_TR), FVAL(KNEW_TR)] = [XOPT + D, F(XOPT + D)], and the model is updated accordingly,
        % but such a model will not be used in the next trust-region iteration, because a geometry step
        % will be invoked to improve the geometry of the interpolation set and update the model again.
        % 6. DELTA has been updated before arriving here: if REDUCE_RHO = FALSE and SHORTD = TRUE, then
        % DELTA was reduced by a factor of 10; if SHORTD = FALSE, then DELTA was updated by TRRAD after
        % the trust-region iteration.
        % 7. If SHORTD = FALSE and KNEW_TR > 0, then XPT has been updated after the trust-region
        % iteration; if RATIO > 0 in addition, then XOPT has been updated as well.

        xdist = sqrt(sum((xpt - repmat(xopt, 1, npt)).^2, 1)');
        bad_trstep = (shortd || ratio < 0.1 || knew_tr == 0);
        improve_geo = (~reduce_rho_1) && (max(xdist) > 2 * delta) && bad_trstep;

        % If all the interpolation points are close to XOPT and the trust region is small, but the
        % trust-region step is "bad" (SHORTD or RATIO <= 0), then we shrink RHO (update the criterion
        % for the "closeness" and SHORTD). REDUCE_RHO_2 corresponds to Box 10 of the NEWUOA paper.
        % N.B.:
        % 1. The definition of REDUCE_RHO_2 is equivalent to the following:
        % REDUCE_RHO_2 = (.NOT. IMPROVE_GEO) .AND. (MAX(DELTA, DNORM) <= RHO) .AND. BAD_TRSTEP
        % 2. The definition of REDUCE_RHO_2 can be moved downward below IF (IMPROVE_GEO) ... END IF.
        % Even though DNORM gets a new value after the geometry step when IMPROVE_GEO = TRUE, this
        % value does not affect REDUCE_RHO_2, because DNORM comes into play only if IMPROVE_GEO = FALSE.
        % 3. DELTA < DNORM may hold due to the update of DELTA.
        bad_trstep = (shortd || ratio <= 0 || knew_tr == 0);
        reduce_rho_2 = (max(xdist) <= 2 * delta) && (max(delta, dnorm) <= rho) && bad_trstep;

        % Comments on BAD_TRSTEP:
        % 1. Powell used different thresholds (<= 0 and < 0.1) for RATIO in the definitions of BAD_TRSTEP
        % above. Unifying them to <= 0 makes little difference to the performance, sometimes worsening,
        % sometimes improving, but never substantially; unifying them to 0.1 makes little difference to
        % the performance.
        % 2. THEORETICALLY, KNEW_TR == 0 implies RATIO < 0, and hence the above definitions of BAD_TRSTEP
        % are mathematically equivalent to (SHORTD .OR. RATIO < 0.1_RP) or (SHORTD .OR. RATIO <= 0).
        % However, KNEW_TR may turn out 0 due to NaN even if RATIO > 0. See SETDROP_TR for details.

        % NEWUOA never sets IMPROVE_GEO and (REDUCE_RHO_1 .OR. REDUCE_RHO_2) to TRUE simultaneously. So
        % the instructions "IF (IMPROVE_GEO) ... END IF" and "IF (REDUCE_RHO_1 .OR. REDUCE_RHO_2)" can
        % be exchanged without changing the algorithm.
        if improve_geo
            % XPT(:, KNEW_GEO) will be dropped (replaced by XOPT + D below).
            [~, knew_geo] = max(xdist);

            % Set DELBAR, which will be used as the trust region radius for the geometry-improving
            % scheme GEOSTEP. We also need it to decide whether to shift XBASE or not.
            % Note that DELTA has been updated before arriving here. See the comments above the
            % definition of IMPROVE_GEO.
            delbar = max(min(0.1 * max(xdist), 0.5 * delta), rho);

            % Shift XBASE if XOPT may be too far from XBASE.
            if delbar * delbar <= 1.0e-3 * norm(xopt)^2
                [bmat, gq, hq, xbase, xopt, xpt] = shiftbase(idz, pq, zmat, bmat, gq, hq, xbase, xopt, xpt, debugflag);
            end

            % Find a step D so that the geometry of XPT will be improved when XPT(:, KNEW_GEO) is
            % replaced by XOPT + D. The GEOSTEP subroutine will call Powell's BIGLAG and BIGDEN. It will
            % also calculate the VLAG and BETA for this D.
            d = geostep(idz, knew_geo, kopt, bmat, delbar, xpt, zmat, debugflag);

            % Calculate VLAG and BETA for D. It makes uses of XOPT, so this is done before updating XOPT.
            [beta, vlag] = vlagbeta(idz, kopt, bmat, d, xpt, zmat, debugflag);

            % Use the current quadratic model to predict the change in F due to the step D.
            qred = calquad(d, gq, hq, pq, xopt, xpt, debugflag);

            % Calculate the next value of the objective function.
            xnew = xopt + d;
            x = xbase + xnew;
            if any(isnan(x))
                f = sum(x);                % Set F to NaN. It is necessary.
                info = infos('nan_x');
                break;
            end
            f = calfun(x);
            nf = nf + 1;
            fmssg(iprint, nf, f, x, solver);
            if maxfhist >= 1
                khist = mod(nf - 1, maxfhist) + 1;
                fhist(khist) = f;
            end
            if maxxhist >= 1
                khist = mod(nf - 1, maxxhist) + 1;
                xhist(:, khist) = x;
            end

            % DNORMSAVE contains the DNORM of the latest 3 function evaluations with the current RHO.
            %------------------------------------------------------------------------------------------%
            % Powell's code does not update DNORM. Therefore, DNORM is the length of last trust-region
            % trial step, which seems inconsistent with what is described in Section 7 (around (7.7)) of
            % the NEWUOA paper. Seemingly we should keep DNORM = ||D|| as we do here. The value of DNORM
            % will be used when defining REDUCE_RHO.
            dnorm = min(delbar, sqrt(norm(d)^2));
            % In theory, DNORM = DELBAR in this case.
            %------------------------------------------------------------------------------------------%
            dnormsav = [dnormsav(2:length(dnormsav)); dnorm];

            % MODERR is the error of the current model in predicting the change in F due to D.
            moderr = f - fopt - qred;
            % MODERRSAVE is the prediction errors of the latest 3 models with the current RHO.
            moderrsav = [moderrsav(2:length(moderrsav)); moderr];

            % Update XOPT and FOPT. Before KOPT is updated, XOPT may differ from XPT(:, KOPT), and FOPT
            % may differ from FVAL(KOPT). Note that the code may exit before KOPT is updated. See below.
            if f < fopt
                xopt = xnew;
                fopt = f;
            end

            % Check whether to exit.
            if isnan(f) || (f == Inf)
                info = infos('nan_inf_f');
                break;
            end
            if f <= ftarget
                info = infos('ftarget_achieved');
                break;
            end
            if nf >= maxfun
                info = infos('maxfun_reached');
                break;
            end

            % Update BMAT, ZMAT and IDZ, so that the KNEW_GEO-th interpolation point can be moved.
            [idz, bmat, zmat] = updateh(knew_geo, beta, vlag, idz, bmat, zmat, debugflag);

            % Update the quadratic model using the updated BMAT, ZMAT, IDZ.
            [gq, hq, pq] = updateq(idz, knew_geo, bmat(:, knew_geo), moderr, zmat, xpt(:, knew_geo), gq, hq, pq, debugflag);

            % Include the new interpolation point.
            xpt(:, knew_geo) = xnew;            % Should be done after UPDATEQ.
            fval(knew_geo) = f;
            if fval(knew_geo) < fval(kopt)
                kopt = knew_geo;
            end
            % KOPT is NOT identical to MINLOC(FVAL). Indeed, if FVAL(KNEW_GEO) = FVAL(KOPT) and
            % KNEW_GEO < KOPT, then MINLOC(FVAL) = KNEW_GEO /= KOPT. Do not change KOPT in this case.
        end        % The procedure of improving geometry ends.

        if reduce_rho_1 || reduce_rho_2
            % The calculations with the current RHO are complete. Pick the next values of RHO and DELTA.
            if rho <= rhoend
                info = infos('small_tr_radius');
                break;
            else
                delta = 0.5 * rho;
                rho_ratio = rho / rhoend;
                if rho_ratio <= 16.0
                    rho = rhoend;
                elseif rho_ratio <= 250.0
                    rho = sqrt(rho_ratio) * rhoend;
                else
                    rho = 0.1 * rho;
                end
                rhomssg(iprint, nf, fopt, rho, xbase + xopt, solver);
                delta = max(delta, rho);
                % DNORMSAVE and MODERRSAVE are corresponding to the latest 3 function evaluations with
                % the current RHO. Update them after reducing RHO.
                dnormsav = Inf(3, 1);
                moderrsav = Inf(length(dnormsav), 1);
            end
        end        % The procedure of reducing RHO ends.

    end    % The iterative procedure ends.

    % Return from the calculation, after another Newton-Raphson step, if it is too short to have been
    % tried before.
    if info == infos('small_tr_radius') && shortd && nf < maxfun
        x = xbase + (xopt + d);
        if any(isnan(x))
            f = sum(x);            % Set F to NaN. It is necessary.
            info = infos('nan_x');
        else
            f = calfun(x);
            nf = nf + 1;
            fmssg(iprint, nf, f, x, solver);
            if maxfhist >= 1
                khist = mod(nf - 1, maxfhist) + 1;
                fhist(khist) = f;
            end
            if maxxhist >= 1
                khist = mod(nf - 1, maxxhist) + 1;
                xhist(:, khist) = x;
            end
        end
    end

    % Choose the [X, F] to return: either the current [X, F] or [XBASE+XOPT, FOPT].
    if isnan(f) || fopt <= f
        x = xbase + xopt;
        f = fopt;
    end

    % Rearrange FHIST and XHIST so that they are in the chronological order.
    if maxfhist >= 1 && maxfhist < nf
        khist = mod(nf - 1, maxfhist) + 1;
        fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)];
    end
    if maxxhist >= 1 && maxxhist < nf
        khist = mod(nf - 1, maxxhist) + 1;
        xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)];
        % N.B.:
        % 1. The result of the array constructor is always a rank-1 array (e.g., vector), no matter what
        % elements are used to construct the array.
        % 2. The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
        % order of Fortran arrays.
        % 3. In MATLAB, `xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)]` does the same thing.
    end

    retmssg(info, iprint, nf, f, x, solver);

end
