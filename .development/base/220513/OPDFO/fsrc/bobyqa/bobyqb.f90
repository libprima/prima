module bobyqb_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the major calculations of BOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Tuesday, July 05, 2022 AM11:26:06
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: bobyqb


contains


subroutine bobyqb(calfun, n, npt, x, xl, xu, rhobeg, rhoend, iprint, &
    & maxfun, xbase, xpt, fval, xopt, gopt, hq, pq, bmat, zmat, ndim, &
    & sl, su, xnew, xalt, d, vlag, w, f, info, ftarget, &
    & nf, xhist, maxxhist, fhist, maxfhist)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist, rangehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
use, non_intrinsic :: info_mod, only : DAMAGING_ROUNDING, SMALL_TR_RADIUS
use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: vm_mod, only : v2m

! Solver-specific modules
use, non_intrinsic :: initialize_mod, only : initialize
use, non_intrinsic :: geometry_mod, only : geostep
use, non_intrinsic :: rescue_mod, only : rescue
use, non_intrinsic :: trustregion_mod, only : trsbox
use, non_intrinsic :: update_mod, only : update

implicit none

! Inputs
procedure(OBJ) :: calfun
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfhist
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: maxxhist
integer(IK), intent(in) :: n
integer(IK), intent(in) :: ndim
integer(IK), intent(in) :: npt
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend
real(RP), intent(in) :: xl(n)
real(RP), intent(in) :: xu(n)

! In-outputs
real(RP), intent(inout) :: sl(n)
real(RP), intent(inout) :: su(n)
real(RP), intent(inout) :: x(n)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: nf
real(RP), intent(out) :: bmat(n, npt + n)
real(RP), intent(out) :: d(n)
real(RP), intent(out) :: f
real(RP), intent(out) :: fhist(maxfhist)
real(RP), intent(out) :: fval(npt)
real(RP), intent(out) :: gopt(n)
real(RP), intent(out) :: hq(n * (n + 1_IK) / 2_IK)
real(RP), intent(out) :: pq(npt)
real(RP), intent(out) :: vlag(npt + n)
real(RP), intent(out) :: w(3_IK * (npt + n))
real(RP), intent(out) :: xalt(n)
real(RP), intent(out) :: xbase(n)
real(RP), intent(out) :: xhist(n, maxxhist)
real(RP), intent(out) :: xnew(n)
real(RP), intent(out) :: xopt(n)
real(RP), intent(out) :: xpt(n, npt)
real(RP), intent(out) :: zmat(npt, npt - n - 1_IK)

! Local variables
real(RP) :: adelt, alpha, bdtest, bdtol, beta, &
&        biglsq, bsumm, cauchy, crvmin, curv, delsq, delta,  &
&        den, denom, densav, diff, diffa, diffb, diffc,     &
&        dist, distsq, dnorm, dsq, dx, errbig, fopt,        &
&        fracsq, frhosq, gisq, gqsq, hdiag,      &
&        pqold, ratio, rho, scaden, summ, summm, summa, summb, summpq,&
&        summw, summz, temp, vquad, xoptsq
integer(IK) :: i, ih, ip, itest, j, jj, jp, k, kbase, knew, &
&           kopt, ksav, nfsav, nh, np, nptm, nresc, ntrits
real(RP) :: bup1(n, npt), bup2(n, n), vtmp(npt), hqm(n, n), gtmp(n)

!     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
!       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT is a TWO-dimensional array that holds the coordinates of the
!       interpolation points relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XOPT is set to the displacement from XBASE of the trust region centre.
!     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
!       this factorization being ZMAT times ZMAT^T, which provides both the
!       correct rank and positive semi-definiTENess.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
!       All the compONEnts of every XOPT are going to satisfy the bounds
!       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
!       XOPT is on a constraint boundary.
!     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
!       vector of variables for the next call of CALFUN. XNEW also satisfies
!       the SL and SU constraints in the way that has just been mentiONEd.
!     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
!       in order to increase the denominator in the updating of UPDATE.
!     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
!     VLAG contains the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     W is a ONE-dimensional array that is used for working space. Its length
!       must be at least 3*NDIM = 3*(NPT+N).
!
!     Set some constants.
!
np = n + 1
nptm = npt - np
nh = (n * np) / 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, with the corresponding values of
!     of NF and KOPT, which are the number of calls of CALFUN so far and the
!     index of the interpolation point at the trust region centre. Then the
!     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
!     less than NPT. GOPT will be updated if KOPT is different from KBASE.
!
call initialize(calfun, n, npt, x, xl, xu, rhobeg, iprint, maxfun, xbase, xpt, &
& fval, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, kopt, f, ftarget, &
& xhist, maxxhist, fhist, maxfhist)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
xopt = xpt(:, kopt)
fopt = fval(kopt)
x = xbase + xopt  ! Set X.
f = fopt  ! Set F.

xoptsq = ZERO
do i = 1, n
    xopt(i) = xpt(i, kopt)
    xoptsq = xoptsq + xopt(i)**2
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: FSAVE is not needed any more. See line number 720.
!      FSAVE=FVAL(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom/Zaikun (on 04-06-2019/07-06-2019):
if (is_nan(f) .or. is_posinf(f)) then
    info = -2
    goto 720
end if
!     By Tom (on 04-06-2019):
!     If F reached the target function, PRELIM will stop and BOBYQB
!     should stop here.
if (f <= ftarget) then
    info = 1
    goto 736
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (nf < npt) then
    info = 3
    goto 720
end if
kbase = 1
!
!     Complete the settings that are required for the iterative procedure.
!
rho = rhobeg
delta = rho
nresc = nf
ntrits = 0
diffa = ZERO
diffb = ZERO
itest = 0
nfsav = nf
!
!     Update GOPT if necessary before the first iteration and after each
!     call of RESCUE that makes a call of CALFUN.
!
20 if (kopt /= kbase) then
    ih = 0
    !do j = 1, n
    !    do i = 1, j
    !        ih = ih + 1
    !        if (i < j) gopt(j) = gopt(j) + hq(ih) * xopt(i)
    !        gopt(i) = gopt(i) + hq(ih) * xopt(j)
    !    end do
    !end do
    hqm = v2m(hq)
    if (nf > npt) then
        !------------------!
        gtmp = gopt; gopt = ZERO
        !------------------!
        do k = 1, npt
            temp = ZERO
            do j = 1, n
                temp = temp + xpt(j, k) * xopt(j)
            end do
            temp = pq(k) * temp
            do i = 1, n
                gopt(i) = gopt(i) + temp * xpt(i, k)
            end do
        end do
        do i = 1, n
            gopt = gopt + xopt(i) * hqm(:, i)
        end do
        !------------------!
        gopt = gtmp + gopt
        !------------------!
    else
        gopt = gopt + matprod(hqm, xopt)
    end if
end if
!
!     Generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     The integer NTRITS is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. If the length
!     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
!     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: For ill-conditiONEd problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRBOX, ALTMOV, or RESCUE is not predictable, and Segmentation Fault or
! infinite cycling may happen. This is because any equality/inequality
! comparison involving NaN returns FALSE, which can lead to uninTENded
! behavior of the code, including uninitialized indices.
!
!   60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,

60 do i = 1, n
    if (gopt(i) /= gopt(i)) then
        info = -3
        goto 720
    end if
end do
do i = 1, nh
    if (hq(i) /= hq(i)) then
        info = -3
        goto 720
    end if
end do
do i = 1, npt
    if (pq(i) /= pq(i)) then
        info = -3
        goto 720
    end if
end do
call trsbox(n, npt, xpt, xopt, gopt, hq, pq, sl, su, delta, xnew, d, &
& w, w(np), w(np + n), w(np + 2 * n), w(np + 3 * n), dsq, crvmin)
dnorm = min(delta, sqrt(dsq))
if (dnorm < HALF * rho) then
    ntrits = -1
    distsq = (TEN * rho)**2
    if (nf <= nfsav + 2) goto 650
!
!     The following choice between labels 650 and 680 depends on whether or
!     not our work with the current RHO seems to be complete. Either RHO is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance HALF*RHO of XOPT.
!
    errbig = max(diffa, diffb, diffc)
    frhosq = 0.125_RP * rho * rho
    if (crvmin > ZERO .and. errbig > frhosq * crvmin) goto 650
    bdtol = errbig / rho
    do j = 1, n
        bdtest = bdtol
        if (xnew(j) == sl(j)) bdtest = w(j)
        if (xnew(j) == su(j)) bdtest = -w(j)
        if (bdtest < bdtol) then
            !---------------------------------!
            ! Zaikun 20220419
            !curv = hq((j + j * j) / 2)
            curv = ZERO
            !---------------------------------!
            do k = 1, npt
                curv = curv + pq(k) * xpt(j, k)**2
            end do
            !---------------------------------!
            ! Zaikun 20220419
            curv = hq((j + j * j) / 2) + curv
            !---------------------------------!
            bdtest = bdtest + HALF * curv * rho
            if (bdtest < bdtol) goto 650
        end if
    end do
    goto 680
end if
ntrits = ntrits + 1
!
!     Severe cancellation is likely to occur if XOPT is too far from XBASE.
!     If the following test holds, then XBASE is shifted so that XOPT becomes
!     ZERO. The appropriate changes are made to BMAT and to the second
!     derivatives of the current model, beginning with the changes to BMAT
!     that do not depend on ZMAT. VLAG is used temporarily for working space.
!
!90 if (dsq <= 1.0E-3_RP * xoptsq) then
90 continue

if (xoptsq >= 1.0E3_RP * dsq) then
    fracsq = 0.25_RP * xoptsq
    summpq = ZERO

    bup2 = ZERO
    do k = 1, npt
        summpq = summpq + pq(k)
        !-----------------------------------------!
        ! Zaikun 20220403
        !summ = -HALF * xoptsq
        summ = ZERO
        do i = 1, n
            summ = summ + xpt(i, k) * xopt(i)
        end do
        !summ = -HALF * xoptsq + summ
        summ = summ - HALF * xoptsq
        !-----------------------------------------!
        w(npt + k) = summ
        temp = fracsq - HALF * summ

        do i = 1, n
            w(i) = bmat(i, k)
            !vlag(i) = summ * xpt(i, k) + temp * xopt(i)
            vlag(i) = summ * (xpt(i, k) - HALF * xopt(i)) + fracsq * xopt(i)

            !ip = npt + i
            do j = 1, i
                !bmat(j, ip) = bmat(j, ip) + w(i) * vlag(j) + vlag(i) * w(j)
                bup2(i, j) = bup2(i, j) + w(i) * vlag(j)
                if (j < i) bup2(j, i) = bup2(j, i) + vlag(i) * w(j)
            end do
        end do
    end do

    !----------------------------------------------------------------!
    ! Zaikun 20220405
    do i = 1, n
        ip = npt + i
        do j = 1, i
            bmat(j, ip) = bmat(j, ip) + (bup2(i, j) + bup2(j, i))
        end do
    end do
    !----------------------------------------------------------------!
!
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
    bup1 = ZERO; bup2 = ZERO
    do jj = 1, nptm
        summz = ZERO
        summw = ZERO
        do k = 1, npt
            summz = summz + zmat(k, jj)
            vlag(k) = w(npt + k) * zmat(k, jj)
            summw = summw + vlag(k)
        end do

        do j = 1, n
            !--------------------------------------------------------------!
            ! Zaikun 20220403
            !summ = (fracsq * summz - HALF * summw) * xopt(j)
            summ = ZERO

            !do k = 1, npt
            !    summ = summ + vlag(k) * xpt(j, k)
            !end do
            !summ = (fracsq * summz - HALF * summw) * xopt(j) + summ
            !--------------------------------------------------------------!
            !--------------------------------------------------------------!
            ! Zaikun 20220405
            do k = 1, npt
                summ = summ + (w(npt + k) * (xpt(j, k) - HALF * xopt(j)) + fracsq * xopt(j)) * zmat(k, jj)
            end do
            !--------------------------------------------------------------!
            w(j) = summ
            do k = 1, npt
                !bmat(j, k) = bmat(j, k) + summ * zmat(k, jj)
                bup1(j, k) = bup1(j, k) + summ * zmat(k, jj)
            end do
        end do
        do i = 1, n
            !ip = i + npt
            temp = w(i)
            do j = 1, i
                !bmat(j, ip) = bmat(j, ip) + temp * w(j)
                bup2(j, i) = bup2(j, i) + temp * w(j)
            end do
        end do
    end do

    !----------------------------------------------------------------!
    ! Zaikun 20220405
    do i = 1, n
        do j = 1, npt
            bmat(i, j) = bmat(i, j) + bup1(i, j)
        end do
    end do

    do i = 1, n
        ip = npt + i
        do j = 1, i
            bmat(j, ip) = bmat(j, ip) + bup2(j, i)
        end do
    end do
    !----------------------------------------------------------------!
!
!     The following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.
!
    ih = 0
    do j = 1, n
        !--------------------------------------------------!
        ! Zaikun 20220403
        !w(j) = -HALF * summpq * xopt(j)
        w(j) = ZERO
        do k = 1, npt
            w(j) = w(j) + pq(k) * xpt(j, k)
            xpt(j, k) = xpt(j, k) - xopt(j)
        end do
        w(j) = -HALF * summpq * xopt(j) + w(j)
        !--------------------------------------------------!
        do i = 1, j
            ih = ih + 1
            !hq(ih) = hq(ih) + w(i) * xopt(j) + xopt(i) * w(j)
            hq(ih) = hq(ih) + (w(i) * xopt(j) + xopt(i) * w(j))
            bmat(j, npt + i) = bmat(i, npt + j)
        end do
    end do
    do i = 1, n
        xbase(i) = xbase(i) + xopt(i)
        xnew(i) = xnew(i) - xopt(i)
        sl(i) = sl(i) - xopt(i)
        su(i) = su(i) - xopt(i)
        xopt(i) = ZERO
    end do
    xoptsq = ZERO
end if
if (ntrits == 0) goto 210
goto 230
!
!     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
!     more expensive than the previous shift, because new matrices BMAT and
!     ZMAT are generated from scratch, which may include the replacement of
!     interpolation points whose positions seem to be causing near linear
!     dependence in the interpolation conditions. Therefore RESCUE is called
!     only if rounding errors have reduced by at least a factor of TWO the
!     denominator of the formula for updating the H matrix. It provides a
!     useful safeguard, but is not invoked in most applications of BOBYQA.
!
190 nfsav = nf
kbase = kopt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29
! See the comments above line number 60.
do i = 1, n
    if (gopt(i) /= gopt(i)) then
        info = -3
        goto 720
    end if
end do
do i = 1, nh
    if (hq(i) /= hq(i)) then
        info = -3
        goto 720
    end if
end do
do i = 1, npt
    if (pq(i) /= pq(i)) then
        info = -3
        goto 720
    end if
end do
do j = 1, n
    do i = 1, ndim
        if (bmat(j, i) /= bmat(j, i)) then
            info = -3
            goto 720
        end if
    end do
end do
do j = 1, nptm
    do i = 1, npt
        if (zmat(i, j) /= zmat(i, j)) then
            info = -3
            goto 720
        end if
    end do
end do

!------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call rescue(calfun, n, npt, xl, xu, iprint, maxfun, xbase, xpt, fval, &
& xopt, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, delta, kopt, &
& vlag, w, w(n + np), w(ndim + np), f, ftarget, &
& xhist, maxxhist, fhist, maxfhist)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------------------------!

!
!     XOPT is updated now in case the branch below to label 720 is taken.
!     Any updating of GOPT occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.
!
xoptsq = ZERO
if (kopt /= kbase) then
    do i = 1, n
        xopt(i) = xpt(i, kopt)
        xoptsq = xoptsq + xopt(i)**2
    end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom/Zaikun (on 04-06-2019/07-06-2019):
if (is_nan(f) .or. is_posinf(f)) then
    info = -2
    goto 720
end if
!     By Tom (on 04-06-2019):
!     If F reached the target function, RESCUE will stop and BOBYQB
!     should stop here.
if (f <= ftarget) then
    info = 1
    goto 736
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if (nf < 0) then
if (nf >= maxfun) then
    nf = maxfun
    info = 3
    goto 720
end if
nresc = nf
if (nfsav < nf) then
    nfsav = nf
    goto 20
end if
if (ntrits > 0) goto 60
!
!     Pick TWO alternative vectors of variables, relative to XBASE, that
!     are suitable as new positions of the KNEW-th interpolation point.
!     Firstly, XNEW is set to the point on a line through XOPT and another
!     interpolation point that minimizes the predicted value of the next
!     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
!     and SU bounds. Secondly, XALT is set to the best feasible point on
!     a constrained version of the Cauchy step of the KNEW-th Lagrange
!     function, the corresponding value of the square of this function
!     being returned in CAUCHY. The choice between these alternatives is
!     going to be made when the denominator is calculated.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Zaikun 23-07-2019:
!  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
!     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
!
!  Although very rare, NaN can sometimes occur in BMAT or ZMAT. If it
!  happens, we terminate the code. See the comments above line number 60.
!  Indeed, if ALTMOV is called with such matrices, then geostep.f90 will
!  encounter a memory error at lines 173--174. This is because the first
!  value of PREDSQ in ALTOMOV (see line 159 of geostep.f90) will be NaN, line
!  164 will not be reached, and hence no value will be assigned to IBDSAV.
!
!  Such an error was observed when BOBYQA was (mistakenly) tested on CUTEst
!  problem CONCON. CONCON is a nonlinearly constrained problem with
!  bounds. By mistake, BOBYQA was called to solve this problem,
!  neglecting all the constraints other than bounds. With only the bound
!  constraints, the objective function turned to be unbounded from
!  below, which led to abnormal values in BMAT (indeed, BETA defined in
!  lines 366--389 took NaN/infinite values).
!
210 do j = 1, n
    do i = 1, ndim
        if (bmat(j, i) /= bmat(j, i)) then
            info = -3
            goto 720
        end if
    end do
end do
do j = 1, nptm
    do i = 1, npt
        if (zmat(i, j) /= zmat(i, j)) then
            info = -3
            goto 720
        end if
    end do
end do
call geostep(n, npt, xpt, xopt, bmat, zmat, ndim, sl, su, kopt, &
& knew, adelt, xnew, xalt, cauchy, w, w(np), w(ndim + 1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
alpha = sum(zmat(knew, :)**2)
do i = 1, n
    d(i) = xnew(i) - xopt(i)
end do
!
!     Calculate VLAG and BETA for the current choice of D. The scalar
!     product of D with XPT(K,.) is going to be held in W(NPT+K) for
!     use when VQUAD is calculated.
!
230 do k = 1, npt
    summa = ZERO
    summb = ZERO
    summ = ZERO
    do j = 1, n
        summa = summa + xpt(j, k) * d(j)
        summb = summb + xpt(j, k) * xopt(j)
        summ = summ + bmat(j, k) * d(j)
    end do
    w(k) = summa * (HALF * summa + summb)
    vlag(k) = summ
    w(npt + k) = summa
end do

!-----------------------------------------!
vtmp = vlag(1:npt)
vlag(1:npt) = ZERO
!-----------------------------------------!
beta = ZERO
do jj = 1, nptm
    summ = ZERO
    do k = 1, npt
        summ = summ + zmat(k, jj) * w(k)
    end do
    beta = beta - summ * summ
    do k = 1, npt
        vlag(k) = vlag(k) + summ * zmat(k, jj)
    end do
end do
!-----------------------------------------!
vlag(1:npt) = vtmp + vlag(1:npt)
!-----------------------------------------!

dsq = ZERO
bsumm = ZERO
dx = ZERO
do j = 1, n
    dsq = dsq + d(j)**2
    summ = ZERO
    do k = 1, npt
        summ = summ + w(k) * bmat(j, k)
    end do
    bsumm = bsumm + summ * d(j)
    jp = npt + j
    do i = 1, n
        summ = summ + bmat(i, jp) * d(i)
    end do
    vlag(jp) = summ
    bsumm = bsumm + summ * d(j)
    dx = dx + d(j) * xopt(j)
end do
!beta = dx * dx + dsq * (xoptsq + dx + dx + HALF * dsq) + beta - bsumm

!bsumm = sum(matprod(bmat(:, npt + 1:npt + n), d) * d(1:n) &
!     & + matprod(bmat(:, 1:npt), w(1:npt)) * d(1:n) &
!     & + matprod(bmat(:, 1:npt), w(1:npt)) * d(1:n))
!beta = dx**2 + dsq * (xoptsq + 2.0_RP * dx + HALF * dsq) + beta - bsumm

!beta = inprod(vlag(1:npt + n), [w(1:npt), d])
!beta = HALF * (inprod(xopt + d, xopt + d)**2 + inprod(xopt, xopt)**2) - inprod(xopt + d, xopt)**2 - beta
beta = dx**2 + dsq * (xoptsq + dx + dx + half * dsq) - inprod(d, vlag(npt + 1:npt + n)) - inprod(w(1:npt), vlag(1:npt))
vlag(kopt) = vlag(kopt) + ONE
!
!     If NTRITS is ZERO, the denominator may be increased by replacing
!     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
!     rounding errors have damaged the chosen denominator.
!
if (ntrits == 0) then
    !denom = vlag(knew)**2 + alpha * beta
    denom = alpha * beta + vlag(knew)**2
    if (denom < cauchy .and. cauchy > ZERO) then
        do i = 1, n
            xnew(i) = xalt(i)
            d(i) = xnew(i) - xopt(i)
        end do
        cauchy = ZERO
        go to 230
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IF (DENOM .LE. HALF*VLAG(KNEW)**2) THEN
    if (.not. (denom > HALF * vlag(knew)**2)) then
!111111111111111111111!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (nf > nresc) goto 190
        !info = 4
        info = DAMAGING_ROUNDING
        goto 720
    end if
!
!     Alternatively, if NTRITS is positive, then set KNEW to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. Again RESCUE may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     KNEW before calculating the next value of the objective function.
!
else
    !delsq = delta * delta
    delsq = delta**2
    scaden = ZERO
    biglsq = ZERO
    knew = 0
    do k = 1, npt
        if (k == kopt) cycle
        hdiag = ZERO
        do jj = 1, nptm
            hdiag = hdiag + zmat(k, jj)**2
        end do
        !den = beta * hdiag + vlag(k)**2
        den = hdiag * beta + vlag(k)**2
        distsq = ZERO
        do j = 1, n
            distsq = distsq + (xpt(j, k) - xopt(j))**2
        end do
        temp = max(ONE, (distsq / delsq)**2)
        if (temp * den > scaden) then
            scaden = temp * den
            knew = k
            denom = den
        end if
        ! -----------------------------------------------------------------------------------------!
        ! Zaikun 20220419
        ! Surprisingly, MEX in MATLAB R2020a with gfortran 9.4.0 behave randomly when evaluating
        ! max(0.0, NaN): sometimes it returns 0.0, sometimes NaN, even if the compilation option is
        ! the same. This makes the output of the next line unpredictable and leads to disagreement
        ! with the modernized version.
        !biglsq = max(biglsq, temp * vlag(k)**2)
        if (temp * vlag(k)**2 > biglsq) biglsq = temp * vlag(k)**2
        ! -----------------------------------------------------------------------------------------!
    end do
    ! -----------------------------------------------------------------------------------------!
    ! Zaikun 20220419
    !if (scaden <= HALF * biglsq) then
    if (.not. scaden > HALF * biglsq) then
        ! -----------------------------------------------------------------------------------------!
        if (nf > nresc) goto 190
        !info = 4
        info = DAMAGING_ROUNDING
        goto 720
    end if
end if
!
!     Put the variables for the next calculation of the objective function
!       in XNEW, with any adjustments for the bounds.
!
!
!     Calculate the value of the objective function at XBASE+XNEW, unless
!       the limit on the number of calculations of F has been reached.
!
360 do i = 1, n
    x(i) = min(max(xl(i), xbase(i) + xnew(i)), xu(i))
    if (xnew(i) == sl(i)) x(i) = xl(i)
    if (xnew(i) == su(i)) x(i) = xu(i)
end do
if (nf >= maxfun) then
    info = 3
    goto 720
end if
nf = nf + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i = 1, n
    if (x(i) /= x(i)) then
        f = x(i) ! set f to nan
        if (nf == 1) then
            fopt = f
            xopt(1:n) = ZERO
        end if
        info = -1
        goto 720
    end if
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------!
call evaluate(calfun, x, f)
call savehist(nf, x, xhist, f, fhist)
!------------------------------------------------------------------------!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom (on 04-06-2019):
if (is_nan(f) .or. is_posinf(f)) then
    if (nf == 1) then
        fopt = f
        xopt(1:n) = ZERO
    end if
    info = -2
    goto 720
end if
!     By Tom (on 04-06-2019):
!     If F achieves the function value, the algorithm exits.
if (f <= ftarget) then
    info = 1
    goto 736
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (ntrits == -1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: FSAVE is not needed any more. See line number 720.
!          FSAVE=F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    goto 720
end if
!
!     Use the quadratic model to predict the change in F due to the step D,
!       and set DIFF to the error of this prediction.
!
fopt = fval(kopt)
vquad = ZERO
ih = 0
do j = 1, n
    temp = ZERO
    do i = 1, n
        if (i <= j) then
            temp = temp + hq(j * (j - 1) / 2 + i) * d(i)
        else
            temp = temp + hq(i * (i - 1) / 2 + j) * d(i)
        end if
    end do
    vquad = vquad + d(j) * (gopt(j) + HALF * temp)

    !vquad = vquad + d(j) * gopt(j)
    !do i = 1, j
    !    ih = ih + 1
    !    temp = d(i) * d(j)
    !    if (i == j) temp = HALF * temp
    !    vquad = vquad + hq(ih) * temp
    !end do
end do
!do k = 1, npt
!    vquad = vquad + HALF * pq(k) * w(npt + k)**2
!end do
vquad = vquad + HALF * inprod(w(npt + 1:2 * npt), pq(1:npt) * w(npt + 1:2 * npt))

diff = f - fopt - vquad
diffc = diffb
diffb = diffa
diffa = abs(diff)
if (dnorm > rho) nfsav = nf
!
!     Pick the next value of DELTA after a trust region step.
!
if (ntrits > 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          IF (VQUAD .GE. ZERO) THEN
    if (.not. (vquad < ZERO)) then
        info = 2
        goto 720
    end if
    ratio = (f - fopt) / vquad
    if (ratio <= TENTH) then
        delta = min(HALF * delta, dnorm)
    else if (ratio <= 0.7_RP) then
        delta = max(HALF * delta, dnorm)
    else
        delta = max(HALF * delta, dnorm + dnorm)
    end if
    if (delta <= 1.5_RP * rho) delta = rho
!
!     Recalculate KNEW and DENOM if the new F is less than FOPT.
!
    if (f < fopt) then
        ksav = knew
        densav = denom
        !delsq = delta * delta
        delsq = delta**2
        scaden = ZERO
        biglsq = ZERO
        knew = 0
        do k = 1, npt
            hdiag = ZERO
            do jj = 1, nptm
                hdiag = hdiag + zmat(k, jj)**2
            end do
            !den = beta * hdiag + vlag(k)**2
            den = hdiag * beta + vlag(k)**2
            distsq = ZERO
            do j = 1, n
                distsq = distsq + (xpt(j, k) - xnew(j))**2
            end do
            temp = max(ONE, (distsq / delsq)**2)
            if (temp * den > scaden) then
                scaden = temp * den
                knew = k
                denom = den
            end if
            ! -----------------------------------------------------------------------------------------!
            ! Zaikun 20220419
            ! Surprisingly, MEX in MATLAB R2020a with gfortran 9.4.0 behave randomly when evaluating
            ! max(0.0, NaN): sometimes it returns 0.0, sometimes NaN, even if the compilation option is
            ! the same. This makes the output of the next line unpredictable and leads to disagreement
            ! between Powell's code and the modernized version.
            !biglsq = max(biglsq, temp * vlag(k)**2)
            if (temp * vlag(k)**2 > biglsq) biglsq = temp * vlag(k)**2
            !------------------------------------------------------------------------------------------!
        end do
        ! -----------------------------------------------------------------------------------------!
        ! Zaikun 20220419
        !if (scaden <= HALF * biglsq) then
        if (.not. scaden > HALF * biglsq) then
            ! -----------------------------------------------------------------------------------------!
            knew = ksav
            denom = densav
        end if
    end if
end if
!
!     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
!     moved. Also update the second derivative terms of the model.
!
call update(n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w)
ih = 0
pqold = pq(knew)
pq(knew) = ZERO
do i = 1, n
    temp = pqold * xpt(i, knew)
    do j = 1, i
        ih = ih + 1
        hq(ih) = hq(ih) + temp * xpt(j, knew)
    end do
end do
do jj = 1, nptm
    temp = diff * zmat(knew, jj)
    do k = 1, npt
        pq(k) = pq(k) + temp * zmat(k, jj)
    end do
end do
!
!     Include the new interpolation point, and make the changes to GOPT at
!     the old XOPT that are caused by the updating of the quadratic model.
!
fval(knew) = f
do i = 1, n
    xpt(i, knew) = xnew(i)
    !w(i) = bmat(i, knew)
    w(i) = ZERO
end do
do k = 1, npt
    summa = ZERO
    do jj = 1, nptm
        !summa = summa + zmat(knew, jj) * zmat(k, jj)
        summa = summa + diff * zmat(knew, jj) * zmat(k, jj)
    end do
    summb = ZERO
    do j = 1, n
        summb = summb + xpt(j, k) * xopt(j)
    end do
    temp = summa * summb
    do i = 1, n
        w(i) = w(i) + temp * xpt(i, k)
    end do
end do
do i = 1, n
    !gopt(i) = gopt(i) + diff * w(i)
    gopt(i) = gopt(i) + diff * bmat(i, knew) + w(i)
end do
!
!     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
!
if (f < fopt) then
    kopt = knew
    xoptsq = ZERO
    ih = 0
    do j = 1, n
        xopt(j) = xnew(j)
        xoptsq = xoptsq + xopt(j)**2
        !do i = 1, j
        !    ih = ih + 1
        !    if (i < j) gopt(j) = gopt(j) + hq(ih) * d(i)
        !    gopt(i) = gopt(i) + hq(ih) * d(j)
        !end do
    end do
    gtmp = gopt; gopt = ZERO
    do k = 1, npt
        temp = ZERO
        do j = 1, n
            temp = temp + xpt(j, k) * d(j)
        end do
        temp = pq(k) * temp
        do i = 1, n
            gopt(i) = gopt(i) + temp * xpt(i, k)
        end do
    end do
    hqm = v2m(hq)
    do i = 1, n
        gopt = gopt + d(i) * hqm(:, i)
    end do
    gopt = gtmp + gopt
end if
!
!     Calculate the parameters of the least Frobenius norm interpolant to
!     the current data, the gradient of this interpolant at XOPT being put
!     into VLAG(NPT+I), I=1,2,...,N.
!
if (ntrits > 0) then
    do k = 1, npt
        vlag(k) = fval(k) - fval(kopt)
        w(k) = ZERO
    end do
    do j = 1, nptm
        summ = ZERO
        do k = 1, npt
            summ = summ + zmat(k, j) * vlag(k)
        end do
        do k = 1, npt
            w(k) = w(k) + summ * zmat(k, j)
        end do
    end do
    do k = 1, npt
        summ = ZERO
        do j = 1, n
            summ = summ + xpt(j, k) * xopt(j)
        end do
        w(k + npt) = w(k)
        w(k) = summ * w(k)
    end do

    gqsq = ZERO
    gisq = ZERO
    do i = 1, n
        summ = ZERO
        do k = 1, npt
            summ = summ + bmat(i, k) * vlag(k)! + xpt(i, k) * w(k)
        end do

        summm = ZERO
        do k = 1, npt
            summm = summm + xpt(i, k) * w(k)
        end do
        summ = summ + summm
        if (xopt(i) == sl(i)) then
            gqsq = gqsq + min(ZERO, gopt(i))**2
            gisq = gisq + min(ZERO, summ)**2
        else if (xopt(i) == su(i)) then
            gqsq = gqsq + max(ZERO, gopt(i))**2
            gisq = gisq + max(ZERO, summ)**2
        else
            gqsq = gqsq + gopt(i)**2
            gisq = gisq + summ**2
        end if
        vlag(npt + i) = summ
    end do
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
    itest = itest + 1
    if (gqsq < TEN * gisq) itest = 0
    if (itest >= 3) then
        do i = 1, max0(npt, nh)
            if (i <= n) gopt(i) = vlag(npt + i)
            if (i <= npt) pq(i) = w(npt + i)
            if (i <= nh) hq(i) = ZERO
            itest = 0
        end do
    end if
end if
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case NTRITS=0 occurs
!     when the new interpolation point was reached by an alternative step.
!
if (ntrits == 0) goto 60
if (f <= fopt + TENTH * vquad) goto 60
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
distsq = max((TWO * delta)**2, (TEN * rho)**2)
650 knew = 0
do k = 1, npt
    summ = ZERO
    do j = 1, n
        summ = summ + (xpt(j, k) - xopt(j))**2
    end do
    if (summ > distsq) then
        knew = k
        distsq = summ
    end if
end do
!
!     If KNEW is positive, then ALTMOV finds alternative new positions for
!     the KNEW-th interpolation point within distance ADELT of XOPT. It is
!     reached via label 90. Otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current RHO are complete.
!
if (knew > 0) then
    dist = sqrt(distsq)
    if (ntrits == -1) then
        delta = min(TENTH * delta, HALF * dist)
        if (delta <= 1.5_RP * rho) delta = rho
    end if
    ntrits = 0
    adelt = max(min(TENTH * dist, delta), rho)
    dsq = adelt * adelt
    goto 90
end if
if (ntrits == -1) goto 680
if (ratio > ZERO) goto 60
if (max(delta, dnorm) > rho) goto 60
!
!     The calculations with the current value of RHO are complete. Pick the
!       next values of RHO and DELTA.
!
680 if (rho > rhoend) then
    delta = HALF * rho
    ratio = rho / rhoend
    if (ratio <= 16.0_RP) then
        rho = rhoend
    else if (ratio <= 250.0_RP) then
        rho = sqrt(ratio) * rhoend
    else
        rho = TENTH * rho
    end if
    delta = max(delta, rho)
    ntrits = 0
    nfsav = nf
    goto 60
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else
    info = SMALL_TR_RADIUS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if
!
!     Return from the calculation, after another Newton-Raphson step, if
!       it is too short to have been tried before.
!
if (ntrits == -1) goto 360
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  720 IF (FVAL(KOPT) .LE. FSAVE) THEN
!  Why update X only when FVAL(KOPT) .LE. FSAVE? This seems INCORRECT,
!  because it may lead to a return with F and X that are not the best
!  available.
720 if (fval(kopt) <= f .or. f /= f) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, n
        x(i) = min(max(xl(i), xbase(i) + xopt(i)), xu(i))
        if (xopt(i) == sl(i)) x(i) = xl(i)
        if (xopt(i) == su(i)) x(i) = xu(i)
    end do
    f = fval(kopt)
end if

736 call rangehist(nf, xhist, fhist)

close (17)
end subroutine bobyqb


end module bobyqb_mod
