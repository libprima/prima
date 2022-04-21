module trustregion_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Thursday, April 21, 2022 PM06:53:00
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trsbox


contains


subroutine trsbox(delta, gopt, hq, pq, sl, su, xopt, xpt, crvmin, d, dsq, gnew, xnew)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, EPS, HALF, DEBUGGING
use, non_intrinsic :: linalg_mod, only : inprod, issymmetric, trueloc
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Inputs
real(RP), intent(in) :: delta
real(RP), intent(in) :: gopt(:)  ! GOPT(N)
real(RP), intent(in) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(in) :: pq(:)  ! PQ(NPT)
real(RP), intent(in) :: sl(:)  ! SL(N)
real(RP), intent(in) :: su(:)  ! SU(N)
real(RP), intent(in) :: xopt(:)  ! XOPT(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! Outputs
real(RP), intent(out) :: crvmin
real(RP), intent(out) :: d(:)  ! D(N)
real(RP), intent(out) :: dsq
real(RP), intent(out) :: gnew(:)  ! GNEW(N)
real(RP), intent(out) :: xnew(:)  ! XNEW(N)

! Local variables
character(len=*), parameter :: srname = 'TRSBOX'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: hred(size(gopt))
real(RP) :: hs(size(gopt))
real(RP) :: s(size(gopt))
real(RP) :: xbdi(size(gopt))
real(RP) :: angbd, angt, beta, bstep, cth, delsq, dhd, dhs,    &
&        dredg, dredsq, ds, ggsav, gredsq,       &
&        qred, rdnext, rdprev, redmax, rednew,       &
&        redsav, resid, sdec, shs, sredg, ssq, stepsq, sth,&
&        stplen, sbound(size(gopt)), temp, tempa, tempb, xsav, xsum(size(xopt)), xtest(size(xopt))
integer(IK) :: i, iact, isav, itcsav, iterc, itermax, iu, &
&           j, k, nact

! Sizes
n = int(size(gopt), kind(n))
npt = int(size(pq), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is n-by-n and symmetric', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(size(xopt) == n, 'SIZE(XOPT) == N', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(size(d) == n, 'SIZE(D) == N', srname)
    call assert(size(gnew) == n, 'SIZE(GNEW) == N', srname)
    call assert(size(xnew) == n, 'SIZE(XNEW) == N', srname)
end if

!
!     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
!       meanings as the corresponding arguments of BOBYQB.
!     DELTA is the trust region radius for the present calculation, which
!       seeks a small value of the quadratic model within distance DELTA of
!       XOPT subject to the bounds on the variables.
!     XNEW will be set to a new vector of variables that is approximately
!       the ONE that minimizes the quadratic model within the trust region
!       subject to the SL and SU constraints on the variables. It satisfies
!       as equations the bounds that become active during the calculation.
!     D is the calculated trial step from XOPT, generated iteratively from an
!       initial value of ZERO. Thus XNEW is XOPT+D after the final iteration.
!     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
!       when D is updated.
!     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
!       set to -1.0, 0.0, or 1.0, the value being nonZERO if and only if the
!       I-th variable has become fixed at a bound, the bound being SL(I) or
!       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
!       information is accumulated during the construction of XNEW.
!     The arrays S, HS and HRED are also used for working space. They hold the
!       current search direction, and the changes in the gradient of Q along S
!       and the reduced D, respectively, where the reduced D is the same as D,
!       except that the components of the fixed variables are ZERO.
!     DSQ will be set to the square of the length of XNEW-XOPT.
!     CRVMIN is set to ZERO if D reaches the trust region boundary. Otherwise
!       it is set to the least curvature of H that occurs in the conjugate
!       gradient searches that are not restricted by any constraints. The
!       value CRVMIN=-ONE is set, however, if all of these searches are
!       constrained.
!
!     A version of the truncated conjugate gradient is applied. If a line
!     search is restricted by a constraint, then the procedure is restarted,
!     the values of the variables that are at their bounds being fixed. If
!     the trust region boundary is reached, then further changes may be made
!     to D, each ONE being in the two dimensional space that is spanned
!     by the current D and the gradient of Q at XOPT+D, staying on the trust
!     region boundary. Termination occurs when the reduction in Q seems to
!     be close to the greatest reduction that can be achieved.
!
!     Set some constants.
!
!
!     The sign of GOPT(I) gives the sign of the change to the I-th variable
!     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
!     or not to fix the I-th variable at ONE of its bounds initially, with
!     NACT being set to the number of fixed variables. D and GNEW are also
!     set for the first iteration. DELSQ is the upper bound on the sum of
!     squares of the free variables. QRED is the reduction in Q so far.
!
iterc = 0
nact = 0
!--------------------------------------------------------------------------------------------------!
iact = 0; itcsav = 0; itermax = n  ! Without this, G95 complains that they are used uninitialized.
!--------------------------------------------------------------------------------------------------!

xbdi = ZERO
xbdi(trueloc(xopt >= su .and. gopt <= 0)) = ONE
xbdi(trueloc(xopt <= sl .and. gopt >= 0)) = -ONE
nact = nact + count(xbdi /= 0)
d = ZERO
gnew = gopt

delsq = delta * delta
qred = ZERO
crvmin = -ONE

! Set the next search direction of the conjugate gradient method. It is the steepest descent
! direction initially and when the iterations are restarted because a variable has just been fixed
! by a bound, and of course the components of the fixed variables are ZERO. ITERMAX is an upper
! bound on the indices of the conjugate gradient iterations.

20 continue
beta = ZERO
30 continue
if (beta == 0) then
    s = -gnew  ! If we are sure that S contain only finite values, we may merge this case into the next.
else
    s = beta * s - gnew
end if
s(trueloc(xbdi /= 0)) = ZERO
stepsq = sum(s**2)
if (stepsq <= ZERO) goto 190

if (beta == ZERO) then
    gredsq = stepsq
    itermax = iterc + n - nact
end if
if (gredsq * delsq <= 1.0E-4_RP * qred * qred) go to 190

! Multiply the search direction by the second derivative matrix of Q and calculate some scalars for
! the choice of steplength. Then set BSTEP to the length of the the step to the trust region
! boundary and STPLEN to the steplength, ignoring the simple bounds.

goto 210

50 continue
resid = delsq - sum(d(trueloc(xbdi == 0))**2)
ds = inprod(d(trueloc(xbdi == 0)), s(trueloc(xbdi == 0)))
shs = inprod(s(trueloc(xbdi == 0)), hs(trueloc(xbdi == 0)))

if (resid <= ZERO) goto 90
temp = sqrt(stepsq * resid + ds * ds)
if (ds < ZERO) then
! Zaikun 20220210: the above line is the original code of Powell. Surprisingly, it works quite
! differently from the following line. Are they different even in precise arithmetic?
! When DS = 0, what should be the simplest (and potentially the stablest) formulation?
! What if we are at the first iteration? BSTEP = DELTA/|D|? See TRSAPP.F90 of NEWUOA.
!if (ds <= ZERO) then  ! Zaikun 20210925
    bstep = (temp - ds) / stepsq
else
    bstep = resid / (temp + ds)
end if
stplen = bstep
if (shs > ZERO) then
    stplen = min(bstep, gredsq / shs)
end if

! Reduce STPLEN if necessary in order to preserve the simple bounds, letting IACT be the index of
! the new constrained variable.
! Zaikun 20220422: Theory and computation differ considerably in the calculation of STPLEN and IACT.
! 1. Theoretically, the WHERE constructs can simplify (S > 0 .and. XTEST > SU) to (S > 0)
! and (S < 0, XTEST < SL) to (S < 0), which will be equivalent to Powell's original code. However,
! overflow will occur due to huge values in SU or SL that indicate the absence of bounds, and
! Fortran compilers will complain. It is not an issue in MATLAB/Python/Julia/R.
! 2. Theoretically, we can also simplify (S > 0 .and. XTEST > SU) to (XTEST > SU). This is because
! the algorithm intends to ensure that SL <= XSUM <= SU, under which the inequality XTEST(I) > SU(I)
! implies S(I) > 0. Numerically, however, XSUM may violate the bounds (slightly) due to rounding
! errors. If we replace (S > 0 .and. XTEST > SU) with (XTEST > SU), then SBOUND(I) will be -Inf when
! SU(I) - XSUM(I) is negative (although tiny) and S(I) is +0 (positively signed zero), which will
! lead to STPLEN = -Inf and IACT = I > 0. This will trigger a restart of the conjugate gradient
! method with DELSQ updated to DELSQ - D(IACT)**2; if D(IACT)**2 << DELSQ, then DELSQ can remain
! unchanged due to rounding, leading to an infinite cycling.
! 3. Theoretically, the WHERE construct corresponding to S > 0 can calculate SBOUND by
! MIN(STPLEN * S, SU - XSUM) / S instead of (SU - XSUM) / S, since this quotient matters only if it
! is less than STPLEN. The motivation is to avoid overflow even without checking XTEST > XU. However,
! such an implementation worsens the performance of BOBYQA significantly in our test on 20220422.
! Why? Note that the conjugate gradient method restarts when IACT > 0. Due to rounding errors,
! MIN(STPLEN * S, SU - XSUM) / S can frequently contain entries less than STPLEN, leading to a
! positive IACT and hence a restart. This turns out harmful to the performance of the algorithm (WHY?).
! It can be rectified in two ways: use MIN(STPLEN, (SU-XSUM) / S) instead of MIN(STPLEN*S, SU-XSUM)/S,
! or set IACT to a positive value only if the minimum of SBOUND is surely less STPLEN, e.g.
! ANY(SBOUND < (ONE-EPS) * STPLEN). The first method does not avoid overflow and makes little sense.
xsum = xopt + d
xtest = xsum + stplen * s
sbound = stplen
where (s > 0 .and. xtest > su)
    sbound = (su - xsum) / s
end where
where (s < 0 .and. xtest < sl)
    sbound = (sl - xsum) / s
end where
!--------------------------------------------------------------------------------------------------!
! The code below is mathematically equivalent to the above but numerically inferior as explained.
!where (s > 0)
!    sbound = min(stplen * s, su - xsum) / s
!end where
!where (s < 0)
!    sbound = max(stplen * s, sl - xsum) / s
!end where
!--------------------------------------------------------------------------------------------------!
sbound(trueloc(is_nan(sbound))) = stplen  ! Needed? No if we are sure that D and S are finite.
iact = 0
if (any(sbound < (ONE - EPS) * stplen)) then
    iact = int(minloc(sbound, dim=1), IK)
    stplen = minval(sbound)
end if
!--------------------------------------------------------------------------------------------------!
! Alternatively, IACT and STPLEN can be calculated as below. We prefer the implementation above:
! 1. The above code is more explicit; in addition, it is more flexible: we can change the condition
! ANY(SBOUND < STPLEN) to ANY(SBOUND < (1 - EPS) * STPLEN) or ANY(SBOUND < (1 + EPS) * STPLEN),
! depending on whether we believe a false positive or a false negative of IACT > 0 is more harmful
! --- according to our test on 20220422, it is the former, as mentioned above.
! 2. The above version is still valid even if we exchange the two lines of IACT and STPLEN.
!iact = int(minloc([stplen, sbound], dim=1), IK) - 1_IK ! This line cannot be exchanged with the next
!stplen = minval([stplen, sbound]) ! This line cannot be exchanged with the last
!--------------------------------------------------------------------------------------------------!

!!MATLAB code for calculating IACT and ALPHA:
!!xsum = xopt + d;
!!sbound = stplen;
!!sbound(s > 0) = (su(s > 0) - xsum(s < 0)) / s(s > 0);
!!sbound(s < 0) = (sl(s < 0) - xsum(s < 0)) / s(s < 0);
!!sbound(isnan(sbound)) = stplen;
!!if any(sbound < stplen)
!!    [stplen, iact] = min(sbound);
!!end

! Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
sdec = ZERO
if (stplen > ZERO) then
    iterc = iterc + 1
    temp = shs / stepsq
    if (iact == 0 .and. temp > ZERO) then
        crvmin = min(crvmin, temp)
        if (crvmin == -ONE) crvmin = temp
    end if
    ggsav = gredsq
    gredsq = ZERO
    do i = 1, n
        gnew(i) = gnew(i) + stplen * hs(i)
        if (xbdi(i) == ZERO) gredsq = gredsq + gnew(i)**2
        d(i) = d(i) + stplen * s(i)
    end do
    sdec = max(stplen * (ggsav - HALF * stplen * shs), ZERO)
    qred = qred + sdec
end if
!----------------------------------------------------------------------!
!if (is_nan(stplen) .or. iterc > 100 * itermax) goto 190 ! Zaikun 20220401
!----------------------------------------------------------------------!
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
if (iact > 0) then
    nact = nact + 1
    call assert(abs(s(iact)) > 0, 'S(IACT) /= 0', srname)
    xbdi(iact) = sign(ONE, s(iact))  !!MATLAB: xbdi(iact) = sign(s(iact))
    !xbdi(iact) = ONE
    !if (s(iact) < ZERO) xbdi(iact) = -ONE
    delsq = delsq - d(iact)**2
    if (delsq <= ZERO) goto 90
    goto 20  ! Zaikun 20220421 Caution: infinite cycling may occur. Fix it!!!
end if
!
!     If STPLEN is less than BSTEP, then either apply another conjugate
!     gradient iteration or RETURN.
!
if (stplen < bstep) then
    if (iterc == itermax) goto 190
    !if (iterc >= itermax) goto 190 ??? Zaikun 20220401
    !----------------------------------------------------------------------------------------------!
    !if (sdec <= 0.01_RP * qred) goto 190  ! An infinite loop to 30 occurred because sdec became NaN
    if (sdec <= 0.01_RP * qred .or. is_nan(sdec) .or. is_nan(qred)) goto 190  ! Zaikun 20220401
    !----------------------------------------------------------------------------------------------!
    beta = gredsq / ggsav
    goto 30  ! Zaikun 20220421 Caution: infinite cycling may occur. Fix it!!!
end if

90 continue

crvmin = ZERO
!
!     Prepare for the alternative iteration by calculating some scalars and
!     by multiplying the reduced D by the second derivative matrix of Q.
!
100 if (nact >= n - 1) goto 190
dredsq = ZERO
dredg = ZERO
gredsq = ZERO
do i = 1, n
    if (xbdi(i) == ZERO) then
        dredsq = dredsq + d(i)**2
        dredg = dredg + d(i) * gnew(i)
        gredsq = gredsq + gnew(i)**2
        s(i) = d(i)
    else
        s(i) = ZERO
    end if
end do
itcsav = iterc
goto 210
!
!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.
!
120 iterc = iterc + 1
!-------------------------------------------------------!
!if (iterc > 100 * itermax) goto 190 ! Zaikun 20220401
!-------------------------------------------------------!
temp = gredsq * dredsq - dredg * dredg
if (temp <= 1.0E-4_RP * qred * qred) goto 190
temp = sqrt(temp)
do i = 1, n
    if (xbdi(i) == ZERO) then
        s(i) = (dredg * d(i) - dredsq * gnew(i)) / temp
    else
        s(i) = ZERO
    end if
end do
! Zaikun 20210926:
!!! Should we calculate S as in TRSAPP of NEWUOA in order to
! make sure that |S| = |D|??? Namely the following:
! S = something, then S = (S/norm(S))*norm(D).
! Also, should exit if the orthogonality of S and D is damaged, or
! S is  not finite.
! See the corresponding part of TRSAPP.
sredg = -temp
!
!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of HALF the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
angbd = ONE
iact = 0
do i = 1, n
    if (xbdi(i) == ZERO) then
        tempa = xopt(i) + d(i) - sl(i)
        tempb = su(i) - xopt(i) - d(i)
        if (tempa <= ZERO) then
            nact = nact + 1
            xbdi(i) = -ONE
            goto 100
        else if (tempb <= ZERO) then
            nact = nact + 1
            xbdi(i) = ONE
            goto 100
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-15: RATIO is never used
!          RATIO=ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ssq = d(i)**2 + s(i)**2
        !------------------------------------------------------------------------------------------!
        !temp = ssq - (xopt(i) - sl(i))**2  ! Overflow can occur due to huge values in SL
        !if (temp > ZERO) then
        !------------------------------------------------------------------------------------------!
        if (xopt(i) - sl(i) < sqrt(ssq)) then
            temp = ssq - (xopt(i) - sl(i))**2
            temp = sqrt(temp) - s(i)
            if (angbd * temp > tempa) then
                angbd = tempa / temp
                iact = i
                xsav = -ONE
            end if
        end if
        !------------------------------------------------------------------------------------------!
        !temp = ssq - (su(i) - xopt(i))**2  ! Overflow can occur due to huge values in SU
        !if (temp > ZERO) then
        !------------------------------------------------------------------------------------------!
        if (su(i) - xopt(i) < sqrt(ssq)) then
            temp = ssq - (su(i) - xopt(i))**2
            temp = sqrt(temp) + s(i)
            if (angbd * temp > tempb) then
                angbd = tempb / temp
                iact = i
                xsav = ONE
            end if
        end if
    end if
end do
!
!     Calculate HHD and some curvatures for the alternative iteration.
!
goto 210
150 shs = ZERO
dhs = ZERO
dhd = ZERO
do i = 1, n
    if (xbdi(i) == ZERO) then
        shs = shs + s(i) * hs(i)
        dhs = dhs + d(i) * hs(i)
        dhd = dhd + d(i) * hred(i)
    end if
end do
!
!     Seek the greatest reduction in Q for a range of equally spaced values
!     of ANGT in [0,ANGBD], where ANGT is the tangent of HALF the angle of
!     the alternative iteration.
!
redmax = ZERO
isav = 0
redsav = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IU=17.0_RP*ANGBD+3.1_RP
iu = int(17.0_RP * angbd + 3.1_RP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i = 1, iu
    angt = angbd * real(i, RP) / real(iu, RP)
    sth = (angt + angt) / (ONE + angt * angt)
    temp = shs + angt * (angt * dhd - dhs - dhs)
    rednew = sth * (angt * dredg - sredg - HALF * sth * temp)
    if (rednew > redmax) then
        redmax = rednew
        isav = i
        rdprev = redsav
    else if (i == isav + 1) then
        rdnext = rednew
    end if
    redsav = rednew
end do
!
!     Return if the reduction is ZERO. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate SDEC.
!
if (isav == 0) goto 190
if (isav < iu) then
    temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext)
    angt = angbd * (real(isav, RP) + HALF * temp) / real(iu, RP)
end if
cth = (ONE - angt * angt) / (ONE + angt * angt)
sth = (angt + angt) / (ONE + angt * angt)
temp = shs + angt * (angt * dhd - dhs - dhs)
sdec = sth * (angt * dredg - sredg - HALF * sth * temp)
if (sdec <= ZERO) goto 190
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
dredg = ZERO
gredsq = ZERO
do i = 1, n
    gnew(i) = gnew(i) + (cth - ONE) * hred(i) + sth * hs(i)
    if (xbdi(i) == ZERO) then
        d(i) = cth * d(i) + sth * s(i)
        dredg = dredg + d(i) * gnew(i)
        gredsq = gredsq + gnew(i)**2
    end if
    hred(i) = cth * hred(i) + sth * hs(i)
end do
qred = qred + sdec
if (iact > 0 .and. isav == iu) then
    nact = nact + 1
    xbdi(iact) = xsav
    goto 100
end if
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
if (sdec > 0.01_RP * qred) goto 120
190 dsq = ZERO
do i = 1, n
    xnew(i) = max(min(xopt(i) + d(i), su(i)), sl(i))
    if (xbdi(i) == -ONE) xnew(i) = sl(i)
    if (xbdi(i) == ONE) xnew(i) = su(i)
    d(i) = xnew(i) - xopt(i)
    dsq = dsq + d(i)**2
end do
return

!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
210 continue
do j = 1, n
    hs(j) = ZERO
    do i = 1, j
        if (i < j) hs(j) = hs(j) + hq(i, j) * s(i)
        hs(i) = hs(i) + hq(i, j) * s(j)
    end do
end do
do k = 1, npt
    if (pq(k) /= ZERO) then
        temp = ZERO
        do j = 1, n
            temp = temp + xpt(j, k) * s(j)
        end do
        temp = temp * pq(k)
        do i = 1, n
            hs(i) = hs(i) + temp * xpt(i, k)
        end do
    end if
end do
if (crvmin /= ZERO) goto 50
if (iterc > itcsav) goto 150
do i = 1, n
    hred(i) = hs(i)
end do
goto 120

end subroutine trsbox


end module trustregion_mod
