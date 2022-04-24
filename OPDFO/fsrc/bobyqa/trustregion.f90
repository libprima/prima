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
! Last Modified: Sunday, April 24, 2022 AM09:44:05
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: trsbox


contains


subroutine trsbox(n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta, &
     &  xnew, d, gnew, xbdi, s, hs, hred, dsq, crvmin)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, HALF
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: vm_mod, only : v2m

implicit none

! Inputs
integer(IK), intent(in) :: n
integer(IK), intent(in) :: npt
real(RP), intent(in) :: delta
real(RP), intent(in) :: gopt(n)
real(RP), intent(in) :: hq(n * (n + 1_IK) / 2_IK)
real(RP), intent(in) :: pq(npt)
real(RP), intent(in) :: sl(n)
real(RP), intent(in) :: su(n)
real(RP), intent(in) :: xopt(n)
real(RP), intent(in) :: xpt(n, npt)

! In-outputs
real(RP), intent(inout) :: crvmin
real(RP), intent(inout) :: d(n)
real(RP), intent(inout) :: dsq
real(RP), intent(inout) :: gnew(n)
real(RP), intent(inout) :: hred(n)
real(RP), intent(inout) :: hs(n)
real(RP), intent(inout) :: s(n)
real(RP), intent(inout) :: xbdi(n)
real(RP), intent(inout) :: xnew(n)

! Local variables
real(RP) :: angbd, angt, beta, bstep, cth, delsq, dhd, dhs,    &
&        dredg, dredsq, ds, ggsav, gredsq,       &
&        qred, rdnext, rdprev, redmax, rednew,       &
&        redsav, resid, sdec, shs, sredg, ssq, stepsq, sth,&
&        stplen, stplensav, temp, tempa, tempb, xsav, xsum(n), sbound(n), hqm(n, n)
integer(IK) :: i, iact, ih, isav, itcsav, iterc, itermax, iu, &
&           j, k, nact
real(RP) :: tang(n), fval(1000)

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
!       except that the compONEnts of the fixed variables are ZERO.
!     DSQ will be set to the square of the length of XNEW-XOPT.
!     CRVMIN is set to ZERO if D reaches the trust region boundary. Otherwise
!       it is set to the least curvature of H that occurs in the conjugate
!       gradient searches that are not restricted by any constraints. The
!       value CRVMIN=-1.0_RP is set, however, if all of these searches are
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-15: SQSTP is never used
!      SQSTP=ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i = 1, n
    xbdi(i) = ZERO
    if (xopt(i) <= sl(i)) then
        if (gopt(i) >= ZERO) xbdi(i) = -ONE
    else if (xopt(i) >= su(i)) then
        if (gopt(i) <= ZERO) xbdi(i) = ONE
    end if
    if (xbdi(i) /= ZERO) nact = nact + 1
    d(i) = ZERO
    gnew(i) = gopt(i)
end do
delsq = delta * delta
qred = ZERO
crvmin = -ONE
!
!     Set the next search direction of the conjugate gradient method. It is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the compONEnts of the fixed variables are ZERO. ITERMAX is an
!     upper bound on the indices of the conjugate gradient iterations.
!
20 beta = ZERO
30 stepsq = ZERO
do i = 1, n
    if (xbdi(i) /= ZERO) then
        s(i) = ZERO
    else if (beta == ZERO) then
        s(i) = -gnew(i)
    else
        s(i) = beta * s(i) - gnew(i)
    end if
    stepsq = stepsq + s(i)**2
end do
if (stepsq == ZERO) goto 190
if (beta == ZERO) then
    gredsq = stepsq
    itermax = iterc + n - nact
end if
if (gredsq * delsq <= 1.0E-4_RP * qred * qred) go to 190
!
!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BSTEP to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.
!
goto 210
50 resid = delsq
ds = ZERO
shs = ZERO
temp = ZERO
do i = 1, n
    if (xbdi(i) == ZERO) then
        !----------------------------------!
        ! Zaikun 20220420
        temp = temp + d(i)**2
        !resid = resid - d(i)**2
        !----------------------------------!
        ds = ds + s(i) * d(i)
        shs = shs + s(i) * hs(i)
    end if
end do
resid = resid - temp

if (resid <= ZERO) goto 90
temp = sqrt(stepsq * resid + ds * ds)
if (ds < ZERO) then
! Zaikun 20220210: the above line is the original code of Powell. Surprisingly, it works quite
! differently from the following line. Are they different even in precise arithmetic?
! When DS = 0, what should be the simplest (and potentially the stablest) formulation?
! What if we are at the first iteration? BSTEP = DELTA/|D|? See TRSAPP.F90 of NEWUOA.
!if (ds <= ZERO) then  ! zaikun 20210925
    bstep = (temp - ds) / stepsq
else
    bstep = resid / (temp + ds)
end if
stplen = bstep
if (shs > ZERO) then
    stplen = min(bstep, gredsq / shs)
end if

!
!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.
!
iact = 0
stplensav = stplen
!sbound = stplen
xsum = xopt + d
do i = 1, n
    if (s(i) /= ZERO) then
        !if (s(i) > ZERO) then
        if (s(i) > ZERO .and. xsum(i) + stplensav * s(i) > su(i)) then
            !temp = min(stplensav * s(i), su(i) - xsum(i)) / s(i)
            temp = (su(i) - xsum(i)) / s(i)
            !sbound(i) = temp
            if (temp < stplen) then
                stplen = temp
                iact = i
            end if
            !else
        else if (s(i) < ZERO .and. xsum(i) + stplensav * s(i) < sl(i)) then
            !temp = max(stplensav * s(i), sl(i) - xsum) / s(i)
            temp = (sl(i) - xsum(i)) / s(i)
            !sbound(i) = temp
            if (temp < stplen) then
                stplen = temp
                iact = i
            end if
        end if
    end if
end do
!
!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
!
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
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
if (iact > 0) then
    nact = nact + 1
    xbdi(iact) = ONE
    if (s(iact) < ZERO) xbdi(iact) = -ONE
    delsq = delsq - d(iact)**2
    if (delsq <= ZERO) goto 90
    goto 20
end if
!
!     If STPLEN is less than BSTEP, then either apply another conjugate
!     gradient iteration or RETURN.
!
if (stplen < bstep) then
    if (iterc == itermax) goto 190
    !----------------------------------------------------------------------------------------------!
    !if (sdec <= 0.01_RP * qred) goto 190
    if (sdec <= 0.01_RP * qred .or. is_nan(sdec) .or. is_nan(qred)) goto 190  ! Zaikun 20220401
    !----------------------------------------------------------------------------------------------!
    beta = gredsq / ggsav
    goto 30
end if
90 crvmin = ZERO
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
! make sure that |S| = |D|??? Namely, do the following:
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
tang = ONE

do i = 1, n
    if (xbdi(i) == ZERO) then
        !tempa = xopt(i) + d(i) - sl(i)
        !tempb = su(i) - xopt(i) - d(i)
        tempa = xopt(i) + d(i) - sl(i)
        tempb = su(i) - (xopt(i) + d(i))
        !if (tempa <= ZERO) then
        if (xopt(i) + d(i) <= sl(i)) then
            nact = nact + 1
            xbdi(i) = -ONE
            goto 100
            !else if (tempb <= ZERO) then
        else if (xopt(i) + d(i) >= su(i)) then
            nact = nact + 1
            xbdi(i) = ONE
            goto 100
        end if
    end if
end do

do i = 1, n
    if (xbdi(i) == ZERO) then
        ssq = d(i)**2 + s(i)**2
        !----------------------------------!
        !tempa = xopt(i) + d(i) - sl(i)
        !tempb = su(i) - xopt(i) - d(i)
        tempa = xopt(i) + d(i) - sl(i)
        tempb = su(i) - (xopt(i) + d(i))
        !----------------------------------!

        !temp = ssq - (xopt(i) - sl(i))**2
        !if (temp > ZERO) then
        !if ((xopt(i) - sl(i))**2 < ssq) then
        if (xopt(i) - sl(i) < sqrt(ssq)) then
            !temp = ssq - (xopt(i) - sl(i))**2
            temp = max(ZERO, ssq - (xopt(i) - sl(i))**2)
            temp = sqrt(temp) - s(i)
            tang(i) = min(tang(i), tempa / temp)
            !if (angbd * temp > tempa) then
            if (angbd > tempa / temp) then
                angbd = tempa / temp
                iact = i
                xsav = -ONE
            end if
        end if

        !temp = ssq - (su(i) - xopt(i))**2
        !if (temp > ZERO) then
        !if ((su(i) - xopt(i))**2 < ssq) then
        if (su(i) - xopt(i) < sqrt(ssq)) then
            !temp = ssq - (su(i) - xopt(i))**2
            temp = max(ZERO, ssq - (su(i) - xopt(i))**2)
            temp = sqrt(temp) + s(i)
            tang(i) = min(tang(i), tempb / temp)
            !if (angbd * temp > tempb) then
            if (angbd > tempb / temp) then
                angbd = tempb / temp
                iact = i
                xsav = ONE
            end if
        end if
    end if
end do
if (any(is_nan(tang))) goto 190
!-------------------------------!
! Zaikun 20220422
!if (angbd <= 0) goto 190
!-------------------------------!
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
iu = int(17.0_RP * angbd + 3.1_RP)
do i = 1, iu
    !angt = angbd * real(i, RP) / real(iu, RP)
    !angt = ZERO + (angbd - ZERO) * real(i, RP) / real(iu, RP)
    angt = ZERO + ((angbd - ZERO) / real(iu, RP)) * real(i, RP)
    if (i == iu) angt = angbd
    sth = (angt + angt) / (ONE + angt * angt)
    temp = shs + angt * (angt * dhd - dhs - dhs)
    rednew = sth * (angt * dredg - sredg - HALF * sth * temp)
    fval(i) = rednew
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
!if (isav == 0) goto 190
if (isav == 0) then
    angt = ZERO
elseif (isav < iu) then
    temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext)
    !angt = angbd * (real(isav, RP) + HALF * temp) / real(iu, RP)
    angt = ZERO + (angbd - ZERO) * (real(isav, RP) + HALF * temp) / real(iu, RP)
else
!-------------------------------!
! Zaikun 20220422
    angt = angbd  ! ISAV = IU
end if
!-------------------------------!
cth = (ONE - angt * angt) / (ONE + angt * angt)
sth = (angt + angt) / (ONE + angt * angt)
temp = shs + angt * (angt * dhd - dhs - dhs)
sdec = sth * (angt * dredg - sredg - HALF * sth * temp)
!if (sdec <= ZERO) goto 190
if (.not. sdec > ZERO) goto 190
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
!if (iact > 0 .and. angt >= angbd) then  ! ??? Changed
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
!--------------------------------------------------------!
! Zaikun 20220422
!ih = 0
!do j = 1, n
!    hs(j) = ZERO
!    do i = 1, j
!        ih = ih + 1
!        if (i < j) hs(j) = hs(j) + hq(ih) * s(i)
!        hs(i) = hs(i) + hq(ih) * s(j)
!    end do
!end do
hs = ZERO
!--------------------------------------------------------!
do k = 1, npt
    !if (pq(k) /= ZERO) then
    temp = ZERO
    do j = 1, n
        temp = temp + xpt(j, k) * s(j)
    end do
    temp = temp * pq(k)
    do i = 1, n
        hs(i) = hs(i) + temp * xpt(i, k)
    end do
    !end if
end do

hqm = v2m(hq)
do i = 1, n
    hs = hs + s(i) * hqm(:, i)
end do

if (crvmin /= ZERO) goto 50
if (iterc > itcsav) goto 150
do i = 1, n
    hred(i) = hs(i)
end do
goto 120

end subroutine trsbox


end module trustregion_mod
