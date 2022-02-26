subroutine lincob(n, npt, m, amat, b, x, rhobeg, rhoend, iprint, maxfun, &
    & xbase, xpt, fval, xsav, xopt, gopt, hq, pq, bmat, zmat, ndim, &
     &  step, sp, xnew, iact, rescon, qfac, rfac, pqw, w, f, info, ftarget)

use, non_intrinsic :: dirty_temporary_mod4powell_mod
implicit real(kind(0.0D0)) (a - h, o - z)
implicit integer(i - n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dimension amat(n, *), b(*), x(*), xbase(*), xpt(npt, *), fval(*), &
     &  xsav(*), xopt(*), gopt(*), hq(*), pq(*), bmat(ndim, *), &
     &  zmat(npt, *), step(*), sp(*), xnew(*), iact(*), rescon(*), &
     &  qfac(n, *), rfac(*), pqw(*), w(*)
!
!     The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
!       identical to the corresponding arguments in SUBROUTINE LINCOA.
!     AMAT is a matrix whose columns are the constraint gradients, scaled
!       so that they have unit length.
!     B contains on entry the right hand sides of the constraints, scaled
!       as above, but later B is modified for variables relative to XBASE.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT contains the interpolation point coordinates relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XSAV holds the best feasible vector of variables so far, without any
!       shift of origin.
!     XOPT is set to XSAV-XBASE, which is the displacement from XBASE of
!       the feasible vector of variables that provides the least calculated
!       F so far, this vector being the current trust region centre.
!     GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of the big inverse matrix H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix
!       of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
!       where the elements of DZ are plus or minus one, as specified by IDZ.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     STEP is employed for trial steps from XOPT. It is also used for working
!       space when XBASE is shifted and in PRELIM.
!     SP is reserved for the scalar products XOPT^T XPT(K,.), K=1,2,...,NPT,
!       followed by STEP^T XPT(K,.), K=1,2,...,NPT.
!     XNEW is the displacement from XBASE of the vector of variables for
!       the current calculation of F, except that SUBROUTINE TRSTEP uses it
!       for working space.
!     IACT is an integer array for the indices of the active constraints.
!     RESCON holds useful information about the constraint residuals. Every
!       nonnegative RESCON(J) is the residual of the J-th constraint at the
!       current trust region centre. Otherwise, if RESCON(J) is negative, the
!       J-th constraint holds as a strict inequality at the trust region
!       centre, its residual being at least |RESCON(J)|; further, the value
!       of |RESCON(J)| is at least the current trust region radius DELTA.
!     QFAC is the orthogonal part of the QR factorization of the matrix of
!       active constraint gradients, these gradients being ordered in
!       accordance with IACT. When NACT is less than N, columns are added
!       to QFAC to complete an N by N orthogonal matrix, which is important
!       for keeping calculated steps sufficiently close to the boundaries
!       of the active constraints.
!     RFAC is the upper triangular part of this QR factorization, beginning
!       with the first diagonal element, followed by the two elements in the
!       upper triangular part of the second column and so on.
!     PQW is used for working space, mainly for storing second derivative
!       coefficients of quadratic functions. Its length is NPT+N.
!     The array W is also used for working space. The required number of
!       elements, namely MAX[M+3*N,2*M+N,2*NPT], is set in LINCOA.
!
!     Set some constants.
!
np = n + 1
nh = (n * np) / 2
nptm = npt - np
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 15-08-2019
! See the comments below line number 210
imprv = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
!       ZMAT and SP for the first iteration. An important feature is that,
!       if the interpolation point XPT(K,.) is not feasible, where K is any
!       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
!       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
!       is set so that XPT(KOPT,.) is the initial trust region centre.
!
call prelim(n, npt, m, amat, b, x, rhobeg, iprint, xbase, xpt, fval, &
& xsav, xopt, gopt, kopt, hq, pq, bmat, zmat, idz, ndim, sp, rescon, &
& step, pqw, w, f, ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom (on 04-06-2019):
if (is_nan(f) .or. is_posinf(f)) then
    fopt = fval(kopt)
    info = -2
    goto 600
end if
!     By Tom/Zaikun (on 04-06-2019/07-06-2019):
!     Note that we should NOT compare F and FTARGET, because X may not
!     be feasible at the exit of PRELIM.
if (fval(kopt) <= ftarget) then
    f = fval(kopt)
    x(1:n) = xsav(1:n)
    info = 1
    goto 616
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Begin the iterative procedure.
!
nf = npt
fopt = fval(kopt)
rho = rhobeg
delta = rho
ifeas = 0
nact = 0
itest = 3
10 knew = 0
nvala = 0
nvalb = 0
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!       to BMAT that do not depend on ZMAT.
!
20 fsave = fopt
xoptsq = zero
do i = 1, n
    xoptsq = xoptsq + xopt(i)**2
end do
if (xoptsq >= 1.0D4 * delta * delta) then
    qoptsq = 0.25D0 * xoptsq
    do k = 1, npt
        sum = zero
        do i = 1, n
            sum = sum + xpt(k, i) * xopt(i)
        end do
        sum = sum - half * xoptsq
        w(npt + k) = sum
        sp(k) = zero
        do i = 1, n
            xpt(k, i) = xpt(k, i) - half * xopt(i)
            step(i) = bmat(k, i)
            w(i) = sum * xpt(k, i) + qoptsq * xopt(i)
            ip = npt + i
            do j = 1, i
                bmat(ip, j) = bmat(ip, j) + step(i) * w(j) + w(i) * step(j)
            end do
        end do
    end do
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
    do k = 1, nptm
        sumz = zero
        do i = 1, npt
            sumz = sumz + zmat(i, k)
            w(i) = w(npt + i) * zmat(i, k)
        end do
        do j = 1, n
            sum = qoptsq * sumz * xopt(j)
            do i = 1, npt
                sum = sum + w(i) * xpt(i, j)
            end do
            step(j) = sum
            if (k < idz) sum = -sum
            do i = 1, npt
                bmat(i, j) = bmat(i, j) + sum * zmat(i, k)
            end do
        end do
        do i = 1, n
            ip = i + npt
            temp = step(i)
            if (k < idz) temp = -temp
            do j = 1, i
                bmat(ip, j) = bmat(ip, j) + temp * step(j)
            end do
        end do
    end do
!
!     Update the right hand sides of the constraints.
!
    if (m > 0) then
        do j = 1, m
            temp = zero
            do i = 1, n
                temp = temp + amat(i, j) * xopt(i)
            end do
            b(j) = b(j) - temp
        end do
    end if
!
!     The following instructions complete the shift of XBASE, including the
!       changes to the parameters of the quadratic model.
!
    ih = 0
    do j = 1, n
        w(j) = zero
        do k = 1, npt
            w(j) = w(j) + pq(k) * xpt(k, j)
            xpt(k, j) = xpt(k, j) - half * xopt(j)
        end do
        do i = 1, j
            ih = ih + 1
            hq(ih) = hq(ih) + w(i) * xopt(j) + xopt(i) * w(j)
            bmat(npt + i, j) = bmat(npt + j, i)
        end do
    end do
    do j = 1, n
        xbase(j) = xbase(j) + xopt(j)
        xopt(j) = zero
        xpt(kopt, j) = zero
    end do
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 21-03-2020
! Exit if BMAT or ZMAT contians NaN
do j = 1, n
    do i = 1, ndim
        if (bmat(i, j) /= bmat(i, j)) then
            info = -3
            goto 600
        end if
    end do
end do
do j = 1, nptm
    do i = 1, npt
        if (zmat(i, j) /= zmat(i, j)) then
            info = -3
            goto 600
        end if
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     In the case KNEW=0, generate the next trust region step by calling
!       TRSTEP, where SNORM is the current trust region radius initially.
!       The final value of SNORM is the length of the calculated step,
!       except that SNORM is zero on return if the projected gradient is
!       unsuitable for starting the conjugate gradient iterations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRSTEP or QMSTEP is not predictable, and Segmentation Fault or
! infinite cycling may happen. This is because any equality/inequality
! comparison involving NaN returns FALSE, which can lead to unintended
! behavior of the code, including uninitialized indices, which can lead
! to segmentation faults.
do j = 1, n
    if (gopt(j) /= gopt(j)) then
        info = -3
        goto 600
    end if
end do
do i = 1, nh
    if (hq(i) /= hq(i)) then
        info = -3
        goto 600
    end if
end do
do i = 1, npt
    if (pq(i) /= pq(i)) then
        info = -3
        goto 600
    end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
delsav = delta
ksave = knew
if (knew == 0) then
    snorm = delta
    do i = 1, n
        xnew(i) = gopt(i)
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call trstep(n, npt, m, amat, xpt, hq, pq, nact, iact, rescon, &
   & qfac, rfac, snorm, step, xnew, w, w(m + 1), pqw, pqw(np), w(m + np))
!
!     A trust region step is applied whenever its length, namely SNORM, is at
!       least HALF*DELTA. It is also applied if its length is at least 0.1999
!       times DELTA and if a line search of TRSTEP has caused a change to the
!       active set. Otherwise there is a branch below to label 530 or 560.
!
    temp = half * delta
    if (xnew(1) >= half) temp = 0.1999D0 * delta
    if (snorm <= temp) then
        delta = half * delta
        if (delta <= 1.4D0 * rho) delta = rho
        nvala = nvala + 1
        nvalb = nvalb + 1
        temp = snorm / rho
        if (delsav > rho) temp = one
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 24-07-2019
!              IF (TEMP .GE. HALF) NVALA=ZERO
!              IF (TEMP .GE. TENTH) NVALB=ZERO
        if (temp >= half) nvala = 0
        if (temp >= tenth) nvalb = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (delsav > rho) goto 530
        if (nvala < 5 .and. nvalb < 3) goto 530
        if (snorm > zero) ksave = -1
        goto 560
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 24-07-2019
!          NVALA=ZERO
!          NVALB=ZERO
    nvala = 0
    nvalb = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Alternatively, KNEW is positive. Then the model step is calculated
!       within a trust region of radius DEL, after setting the gradient at
!       XBASE and the second derivative parameters of the KNEW-th Lagrange
!       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
!
else
    del = dmax1(tenth * delta, rho)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: See the comments below line number 140
!          DO 160 I=1,N
!  160     W(I)=BMAT(KNEW,I)
    do i = 1, n
        w(i) = bmat(knew, i)
        if (w(i) /= w(i)) then
            info = -3
            goto 600
        end if
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1, npt
        pqw(k) = zero
    end do
    do j = 1, nptm
        temp = zmat(knew, j)
        if (j < idz) temp = -temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: See the comments below line number 140
! Note that the data in PQW is used in QMSTEP below
!          DO 180 K=1,NPT
!  180     PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
        do k = 1, npt
            pqw(k) = pqw(k) + temp * zmat(k, j)
            if (pqw(k) /= pqw(k)) then
                info = -3
                goto 600
            end if
        end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: B is never used in QMSTEP
!          CALL QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,RESCON,
    call qmstep(n, npt, m, amat, xpt, xopt, nact, iact, rescon, &
&   qfac, kopt, knew, del, step, w, pqw, w(np), w(np + m), ifeas)
end if
!
!     Set VQUAD to the change to the quadratic model when the move STEP is
!       made from XOPT. If STEP is a trust region step, then VQUAD should be
!       negative. If it is nonnegative due to rounding errors in this case,
!       there is a branch to label 530 to try to improve the model.
!
vquad = zero
ih = 0
do j = 1, n
    vquad = vquad + step(j) * gopt(j)
    do i = 1, j
        ih = ih + 1
        temp = step(i) * step(j)
        if (i == j) temp = half * temp
        vquad = vquad + temp * hq(ih)
    end do
end do
do k = 1, npt
    temp = zero
    do j = 1, n
        temp = temp + xpt(k, j) * step(j)
        sp(npt + k) = temp
    end do
    vquad = vquad + half * pq(k) * temp * temp
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 15-08-2019
! Although very rarely, with the original code, an infinite loop can occur
! in the following scenario.
! Suppose that, at an certain iteration,
! KNEW = 0, SNORM > 0.5*DELTA > RHO, VQUAD >= 0, and
! sum_{K=1}^NPT ||XPT(K,:)-XOPT(:)||^2 < DELTA^2
! (i.e., DELTA is large and SNORM is not small, yet VQUAD >= 0 due to
! rounding errors and XPT are not far from XOPT).
! Then the program will goto 530 and then goto 20, where XBASE may be
! shifted to the current best point, in the hope of reducing rounding
! errors and 'improve' the model. Afterwards, another trust region step
! is produced by the 'improved' model. Note that DELTA remains unchanged
! in this process. If the new trust region step turns out to satisfy
! SNORM > 0.5*DELTA and VQUAD >= 0 again (i.e., the 'improved' model
! still suffers from rounding errors), then the program will goto 530
! and then goto 20, where shifting will not happen because either XBASE
! was already shifted to the current best point in last step, or XBASE
! is close to the current best point. Consequently, the model will
! remain unchanged, and produce the same trust region step again. This
! leads to an infinite loop.
! The infinite loop did happen when the MATLAB interface was applied to
! min atan(x+100) s.t. x<=-99 (x0=-99, npt=3, rhobeg=1, rhoend=1e-6).
! The problem does not exist in NEWUOA or BOBYQA, where the program will
! exit immediately when VQUAD >= 0.
! To prevent such a loop, here we use IMPRV to record whether the path
! 530 --> 20 has already happened for last trust region step. IMPRV=1
! implies that last trust region step satisfies VQUAD >= 0 and followed
! 530 --> 20. With IMPRV=1, if VQUAD is again nonnegative for the new trust
! region step, we should not goto 530 but goto 560, where IMPRV will be
! set to 0 and DELTA will be reduced. Otherwise, an infinite loop would happen.
!      IF (KSAVE .EQ. 0 .AND. VQUAD .GE. ZERO) GOTO 530
if (ksave == 0 .and. .not. (vquad < zero)) then
    if (imprv == 1) then
        goto 560
    else
        imprv = 1
        goto 530
    end if
else
    imprv = 0
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the next value of the objective function. The difference
!       between the actual new value of F and the value predicted by the
!       model is recorded in DIFF.
!
220 nf = nf + 1
if (nf > maxfun) then
    nf = nf - 1
    info = 3
    goto 600
end if
xdiff = zero
do i = 1, n
    xnew(i) = xopt(i) + step(i)
    x(i) = xbase(i) + xnew(i)
    xdiff = xdiff + (x(i) - xsav(i))**2
end do
xdiff = dsqrt(xdiff)
if (ksave == -1) xdiff = rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IF (XDIFF .LE. TENTH*RHO .OR. XDIFF .GE. DELTA+DELTA) THEN
if (.not. (xdiff > tenth * rho .and. xdiff < delta + delta)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ifeas = 0
    info = 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    goto 600
end if
if (ksave <= 0) ifeas = 1
f = dfloat(ifeas)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i = 1, n
    if (x(i) /= x(i)) then
        f = x(i) ! set f to nan
        if (nf == 1) then
            fopt = f
            xopt(1:n) = zero
        end if
        info = -1
        goto 600
    end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call calfun(n, x, f)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom (on 04-06-2019):
if (is_nan(f) .or. is_posinf(f)) then
    if (nf == 1) then
        fopt = f
        xopt(1:n) = zero
    end if
    info = -2
    goto 600
end if
if (ksave == -1) then
    info = 0
    goto 600
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
diff = f - fopt - vquad
!
!     If X is feasible, then set DFFALT to the difference between the new
!       value of F and the value predicted by the alternative model.
!
if (ifeas == 1 .and. itest < 3) then
    do k = 1, npt
        pqw(k) = zero
        w(k) = fval(k) - fval(kopt)
    end do
    do j = 1, nptm
        sum = zero
        do i = 1, npt
            sum = sum + w(i) * zmat(i, j)
        end do
        if (j < idz) sum = -sum
        do k = 1, npt
            pqw(k) = pqw(k) + sum * zmat(k, j)
        end do
    end do
    vqalt = zero
    do k = 1, npt
        sum = zero
        do j = 1, n
            sum = sum + bmat(k, j) * step(j)
        end do
        vqalt = vqalt + sum * w(k)
        vqalt = vqalt + pqw(k) * sp(npt + k) * (half * sp(npt + k) + sp(k))
    end do
    dffalt = f - fopt - vqalt
end if
if (itest == 3) then
    dffalt = diff
    itest = 0
end if
!
!     Pick the next value of DELTA after a trust region step.
!
if (ksave == 0) then
    ratio = (f - fopt) / vquad
    if (ratio <= tenth) then
        delta = half * delta
    else if (ratio <= 0.7D0) then
        delta = dmax1(half * delta, snorm)
    else
        temp = dsqrt(2.0D0) * delta
        delta = dmax1(half * delta, snorm + snorm)
        delta = dmin1(delta, temp)
    end if
    if (delta <= 1.4D0 * rho) delta = rho
end if
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!       can be moved. If STEP is a trust region step, then KNEW is zero at
!       present, but a positive value is picked by subroutine UPDATE.
!
call update(n, npt, xpt, bmat, zmat, idz, ndim, sp, step, kopt, knew, pqw, w)
if (knew == 0) then
    info = 9
    goto 600
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 19-03-2020
! Exit if BMAT or ZMAT contians NaN
do j = 1, n
    do i = 1, ndim
        if (bmat(i, j) /= bmat(i, j)) then
            info = -3
            goto 600
        end if
    end do
end do
do j = 1, nptm
    do i = 1, npt
        if (zmat(i, j) /= zmat(i, j)) then
            info = -3
            goto 600
        end if
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     If ITEST is increased to 3, then the next quadratic model is the
!       one whose second derivative matrix is least subject to the new
!       interpolation conditions. Otherwise the new model is constructed
!       by the symmetric Broyden method in the usual way.
!
if (ifeas == 1) then
    itest = itest + 1
    if (dabs(dffalt) >= tenth * dabs(diff)) itest = 0
end if
!
!     Update the second derivatives of the model by the symmetric Broyden
!       method, using PQW for the second derivative parameters of the new
!       KNEW-th Lagrange function. The contribution from the old parameter
!       PQ(KNEW) is included in the second derivative matrix HQ. W is used
!       later for the gradient of the new KNEW-th Lagrange function.
!
if (itest < 3) then
    do k = 1, npt
        pqw(k) = zero
    end do
    do j = 1, nptm
        temp = zmat(knew, j)
        if (temp /= zero) then
            if (j < idz) temp = -temp
            do k = 1, npt
                pqw(k) = pqw(k) + temp * zmat(k, j)
            end do
        end if
    end do
    ih = 0
    do i = 1, n
        w(i) = bmat(knew, i)
        temp = pq(knew) * xpt(knew, i)
        do j = 1, i
            ih = ih + 1
            hq(ih) = hq(ih) + temp * xpt(knew, j)
        end do
    end do
    pq(knew) = zero
    do k = 1, npt
        pq(k) = pq(k) + diff * pqw(k)
    end do
end if
!
!     Include the new interpolation point with the corresponding updates of
!       SP. Also make the changes of the symmetric Broyden method to GOPT at
!       the old XOPT if ITEST is less than 3.
!
fval(knew) = f
sp(knew) = sp(kopt) + sp(npt + kopt)
ssq = zero
do i = 1, n
    xpt(knew, i) = xnew(i)
    ssq = ssq + step(i)**2
end do
sp(npt + knew) = sp(npt + kopt) + ssq
if (itest < 3) then
    do k = 1, npt
        temp = pqw(k) * sp(k)
        do i = 1, n
            w(i) = w(i) + temp * xpt(k, i)
        end do
    end do
    do i = 1, n
        gopt(i) = gopt(i) + diff * w(i)
    end do
end if
!
!     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
!       least calculated value so far with a feasible vector of variables.
!
if (f < fopt .and. ifeas == 1) then
    fopt = f
    do j = 1, n
        xsav(j) = x(j)
        xopt(j) = xnew(j)
    end do
    kopt = knew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         By Tom (on 04-06-2019):
    if (fopt <= ftarget) then
        info = 1
        goto 616
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    snorm = dsqrt(ssq)
    do j = 1, m
        if (rescon(j) >= delta + snorm) then
            rescon(j) = snorm - rescon(j)
        else
            rescon(j) = rescon(j) + snorm
            if (rescon(j) + delta > zero) then
                temp = b(j)
                do i = 1, n
                    temp = temp - xopt(i) * amat(i, j)
                end do
                temp = dmax1(temp, zero)
                if (temp >= delta) temp = -temp
                rescon(j) = temp
            end if
        end if
    end do
    do k = 1, npt
        sp(k) = sp(k) + sp(npt + k)
    end do
!
!     Also revise GOPT when symmetric Broyden updating is applied.
!
    if (itest < 3) then
        ih = 0
        do j = 1, n
            do i = 1, j
                ih = ih + 1
                if (i < j) gopt(j) = gopt(j) + hq(ih) * step(i)
                gopt(i) = gopt(i) + hq(ih) * step(j)
            end do
        end do
        do k = 1, npt
            temp = pq(k) * sp(npt + k)
            do i = 1, n
                gopt(i) = gopt(i) + temp * xpt(k, i)
            end do
        end do
    end if
end if
!
!     Replace the current model by the least Frobenius norm interpolant if
!       this interpolant gives substantial reductions in the predictions
!       of values of F at feasible points.
!
if (itest == 3) then
    do k = 1, npt
        pq(k) = zero
        w(k) = fval(k) - fval(kopt)
    end do
    do j = 1, nptm
        sum = zero
        do i = 1, npt
            sum = sum + w(i) * zmat(i, j)
        end do
        if (j < idz) sum = -sum
        do k = 1, npt
            pq(k) = pq(k) + sum * zmat(k, j)
        end do
    end do
    do j = 1, n
        gopt(j) = zero
        do i = 1, npt
            gopt(j) = gopt(j) + w(i) * bmat(i, j)
        end do
    end do
    do k = 1, npt
        temp = pq(k) * sp(k)
        do i = 1, n
            gopt(i) = gopt(i) + temp * xpt(k, i)
        end do
    end do
    do ih = 1, nh
        hq(ih) = zero
    end do
end if
!
!     If a trust region step has provided a sufficient decrease in F, then
!       branch for another trust region calculation. Every iteration that
!       takes a model step is followed by an attempt to take a trust region
!       step.
!
knew = 0
if (ksave > 0) goto 20
if (ratio >= tenth) goto 20
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
530 distsq = dmax1(delta * delta, 4.0D0 * rho * rho)
do k = 1, npt
    sum = zero
    do j = 1, n
        sum = sum + (xpt(k, j) - xopt(j))**2
    end do
    if (sum > distsq) then
        knew = k
        distsq = sum
    end if
end do
!
!     If KNEW is positive, then branch back for the next iteration, which
!       will generate a "model step". Otherwise, if the current iteration
!       has reduced F, or if DELTA was above its lower bound when the last
!       trust region step was calculated, then try a "trust region" step
!       instead.
!
if (knew > 0) goto 20
knew = 0
if (fopt < fsave) goto 20
if (delsav > rho) goto 20
!
!     The calculations with the current value of RHO are complete.
!       Pick the next value of RHO.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 15-08-2019
! See the comments below line number 210
!  560 IF (RHO .GT. RHOEND) THEN
560 imprv = 0
if (rho > rhoend) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    delta = half * rho
    if (rho > 250.0D0 * rhoend) then
        rho = tenth * rho
    else if (rho <= 16.0D0 * rhoend) then
        rho = rhoend
    else
        rho = dsqrt(rho * rhoend)
    end if
    delta = dmax1(delta, rho)
    goto 10
end if
!
!     Return from the calculation, after branching to label 220 for another
!       Newton-Raphson step if it has not been tried before.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (ksave == -1) goto 220
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  600 IF (FOPT .LE. F .OR. IFEAS .EQ. 0) THEN
600 if (fopt <= f .or. ifeas == 0 .or. f /= f) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, n
        x(i) = xsav(i)
    end do
    f = fopt
end if
616 w(1) = f
w(2) = dfloat(nf) + half
return
end
