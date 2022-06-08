module initialize_mod
!--------------------------------------------------------------------------------------------------!
! This module performs the initialization of UOBYQA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
!
! Started: July 2020
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Last Modified: Wednesday, June 08, 2022 PM11:01:49
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initalize


contains


subroutine initalize(calfun, iprint, maxfun, ftarget, rhobeg, x0, kopt, nf, fhist, fopt, xbase, &
    & xhist, xpt, pq, pl, info)

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan, is_posinf
use, non_intrinsic :: info_mod, only : INFO_DFT
use, non_intrinsic :: linalg_mod, only : eye, trueloc, linspace
use, non_intrinsic :: output_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: info_mod, only : NAN_INF_X, NAN_INF_F, NAN_MODEL, FTARGET_ACHIEVED, MAXFUN_REACHED

implicit none

! Inputs
procedure(OBJ) :: calfun  ! N.B.: INTENT cannot be specified if a dummy procedure is not a POINTER
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: x0(:)   ! X0(N)

! Outputs
integer(IK), intent(out) :: info
integer(IK), intent(out) :: kopt
integer(IK), intent(out) :: nf
real(RP), intent(out) :: fhist(:)   ! FHIST(MAXFHIST)
real(RP), intent(out) :: fopt
real(RP), intent(out) :: xbase(:)   ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)    ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(out) :: pl(:, :)  ! PL((N + 1) * (N + 2) / 2, (N + 1) * (N + 2) / 2 - 1)
real(RP), intent(out) :: pq(:)  ! PQ((N + 1) * (N + 2) / 2 - 1)

! Local variables
character(len=*), parameter :: solver = 'UOBYQA'
character(len=*), parameter :: srname = 'INITXF'
integer(IK) :: k
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
integer(IK) :: subinfo
logical :: evaluated(size(xpt, 2))
real(RP) :: f, fval(size(xpt, 2))
real(RP) :: x(size(x0))

integer(IK) :: iw, ih, ip, iq, jj, kk, k0(size(x0)), k1(size(x0))
real(RP) :: rho, rhosq, fbase, xw(size(x0)), d(size(x0)), temp, tempa

n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = max(maxxhist, maxfhist)

rho = rhobeg
rhosq = rho**2

!====================!
! Calculation starts !
!====================!

! Initialize INFO to the default value. At return, an INFO different from this value will indicate
! an abnormal return.
info = INFO_DFT

! Initialize XBASE to X0.
xbase = x0

! EVALUATED is a boolean array with EVALUATED(I) indicating whether the function value of the I-th
! interpolation point has been evaluated. We need it for a portable counting of the number of
! function evaluations, especially if the loop is conducted asynchronously. However, the loop here
! is not fully parallelizable if NPT>2N+1, as the definition XPT(;, 2N+2:end) involves FVAL(1:2N+1).
evaluated = .false.
! Initialize FVAL to HUGENUM. Otherwise, compilers may complain that FVAL is not (completely)
! initialized if the initialization aborts due to abnormality (see CHECKEXIT).
fval = HUGENUM

! Initialization. NF is the number of function calculations so far. The least function value so far,
! the corresponding X and its index are noted in FOPT, XOPT, and KOPT respectively.
info = INFO_DFT
xpt = ZERO
pl = ZERO

!k0 = linspace(2_IK, n + 1_IK, n); k1 = k0 + n
k0 = linspace(2_IK, 2_IK * n, n); k1 = k0 + 1_IK

do k = 1, 1
    x = xpt(:, k) + xbase
    call evaluate(calfun, x, f)
    evaluated(k) = .true.
    fval(k) = f
    call fmsg(solver, iprint, k, f, x)
    ! Save X and F into the history.
    call savehist(k, x, xhist, f, fhist)
    ! Check whether to exit.
    subinfo = checkexit(maxfun, k, f, ftarget, x)
    if (subinfo /= INFO_DFT) then
        info = subinfo
        exit
    end if
end do

if (info == INFO_DFT) then
    do jj = 1, n
        !kk = 2_IK * jj
        k = k0(jj)
        xpt(jj, k) = rho
        x = xpt(:, k) + xbase
        call evaluate(calfun, x, f)
        evaluated(k) = .true.
        fval(k) = f
        call fmsg(solver, iprint, k, f, x)
        ! Save X and F into the history.
        call savehist(k, x, xhist, f, fhist)
        ! Check whether to exit.
        subinfo = checkexit(maxfun, k, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        k = k1(jj)
        if (fval(k0(jj)) < fval(1)) then
            xpt(jj, k) = TWO * rho
        else
            xpt(jj, k) = -rho
        end if
        x = xpt(:, k) + xbase
        call evaluate(calfun, x, f)
        evaluated(k) = .true.
        fval(k) = f
        call fmsg(solver, iprint, k, f, x)
        ! Save X and F into the history.
        call savehist(k, x, xhist, f, fhist)
        ! Check whether to exit.
        subinfo = checkexit(maxfun, k, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end do
end if

xw = -rho
xw(trueloc(fval(k0) < fval(1))) = rho

if (info == INFO_DFT) then
    ip = 0
    iq = 2
    do k = 2_IK * n + 2_IK, npt
        ! Pick the shift from XBASE to the next initial interpolation point that provides
        ! off-diagonal second derivatives.
        ip = ip + 1
        if (ip == iq) then
            ip = 1
            iq = iq + 1
        end if
        ! N.B.: XPT(2, K+1) is accessed by XPT(IQ, K+1) even if N = 1.
        xpt([ip, iq], k) = xw([ip, iq])
        x = xpt(:, k) + xbase
        call evaluate(calfun, x, f)
        evaluated(k) = .true.
        fval(k) = f
        call fmsg(solver, iprint, k, f, x)
        ! Save X and F into the history.
        call savehist(k, x, xhist, f, fhist)
        ! Check whether to exit.
        subinfo = checkexit(maxfun, k, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if
    end do
end if

nf = int(count(evaluated), kind(nf))  !!MATLAB: nf = sum(evaluated);
kopt = int(minloc(fval, mask=evaluated, dim=1), kind(kopt))
fopt = fval(kopt)
!!MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first')


if (all(evaluated)) then
    fbase = fval(1)
    do jj = 1, n
        kk = k0(jj); k = k1(jj)
        ! Form the gradient and diagonal second derivatives of the quadratic model and Lagrange functions.
        ih = n + jj * (jj + 1_IK) / 2_IK
        if (xpt(jj, k) > 0) then  ! XPT(JJ, K) = 2*RHO
            d(jj) = (fbase + fval(k) - TWO * fval(kk)) / rhosq
            pq(jj) = (4.0_RP * fval(kk) - 3.0_RP * fbase - fval(k)) / (TWO * rho)
        else  ! XPT(JJ, K) = -RHO
            d(jj) = (fval(kk) + fval(k) - TWO * fbase) / rhosq
            pq(jj) = (fval(kk) - fval(k)) / (TWO * rho)
        end if
        pq(ih) = d(jj)
    end do

! Form the off-diagonal second derivatives of the initial quadratic model.
    ih = n + 1
    ip = 0
    iq = 2
    do k = 2_IK * n + 2_IK, npt
        ip = ip + 1
        if (ip == iq) then
            ih = ih + 1
            ip = 1
            iq = iq + 1
        end if
        ih = ih + 1
        temp = ONE / (xw(ip) * xw(iq))
        tempa = fval(k) - fbase - xw(ip) * pq(ip) - xw(iq) * pq(iq)
        ! N.B.: D(2) is accessed by D(IQ) even if N = 1.
        pq(ih) = (tempa - HALF * rhosq * (d(ip) + d(iq))) * temp
    end do
end if


if (all(evaluated)) then
    do jj = 1, n
        kk = k0(jj); k = k1(jj)
        ! Form the gradient and diagonal second derivatives of the quadratic model and Lagrange functions.
        ih = n + jj * (jj + 1_IK) / 2_IK
        if (xpt(jj, k) > 0) then  ! XPT(JJ, K) = 2*RHO
            pl(1, jj) = -1.5_RP / rho
            pl(1, ih) = ONE / rhosq
            pl(kk, jj) = TWO / rho  ! Should be moved out of the loop
            pl(kk, ih) = -TWO / rhosq  ! Should be moved out of the loop
        else  ! XPT(JJ, K) = -RHO
            pl(1, ih) = -TWO / rhosq
            pl(kk, jj) = HALF / rho  ! Should be moved out of the loop
            pl(kk, ih) = ONE / rhosq  ! Should be moved out of the loop
        end if
        pl(k, jj) = -HALF / rho
        pl(k, ih) = ONE / rhosq
    end do

! Form the off-diagonal second derivatives of the initial quadratic model.
    ih = n + 1
    ip = 0
    iq = 2
    do k = 2_IK * n + 2_IK, npt
        ip = ip + 1
        if (ip == iq) then
            ih = ih + 1
            ip = 1
            iq = iq + 1
        end if
        ih = ih + 1
        temp = ONE / (xw(ip) * xw(iq))
        pl(1, ih) = temp
        pl(k, ih) = temp

        if (xw(ip) < 0) then
            iw = 2_IK * ip + 1_IK
        else
            iw = 2_IK * ip
        end if
        pl(iw, ih) = -temp

        if (xw(iq) < 0) then
            iw = 2_IK * iq + 1_IK
        else
            iw = 2_IK * iq
        end if
        pl(iw, ih) = -temp
    end do
end if

end subroutine

end module initialize_mod
