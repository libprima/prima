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
! Last Modified: Wednesday, June 08, 2022 PM01:57:21
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
use, non_intrinsic :: linalg_mod, only : eye
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

integer(IK) :: iw, ih, ip, iq, i, j, kk
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
j = 0
ih = n
! In the following loop, FPLUS is set to F(X + RHO*e_I) when NF = 2*I, and the value of FPLUS is
! used subsequently when NF = 2*I + 1.
!fval(kk) = ZERO ! This initial value is not used but to entertain the Fortran compilers.
do k = 1, 2_IK * n + 1_IK
    kk = 2 * (k / 2)
    ! Pick the shift from XBASE to the next initial interpolation point that provides diagonal
    ! second derivatives.
    if (k > 1) then
        if (modulo(k, 2_IK) == 1_IK) then
            if (fval(kk) < fbase) then
                xw(j) = rho
                xpt(j, k) = TWO * rho
            else
                xw(j) = -rho
                xpt(j, k) = -rho
            end if
        else
            j = k / 2_IK
            xpt(j, k) = rho
        end if
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

    if (k == 1) then
        fbase = f
    end if

    if (modulo(k, 2_IK) == 0_IK) then
        fval(kk) = f  ! FPLUS = F(X + RHO*e_I) with I = NF/2.
        cycle
    end if
end do

do k = 1, 2_IK * n + 1_IK
    j = k / 2_IK
    ! Form the gradient and diagonal second derivatives of the quadratic model and Lagrange functions.
    if (j >= 1 .and. k >= 3) then  ! NF >= 3 is implied by J >= 1. We prefer to impose it explicitly.
        ih = ih + j
        if (xpt(j, k) > 0) then  ! XPT(J, NF) = 2*RHO
            pq(j) = (4.0_RP * fval(kk) - 3.0_RP * fbase - fval(k)) / (TWO * rho)
            d(j) = (fbase + fval(k) - TWO * fval(kk)) / rhosq
            pl(1, j) = -1.5_RP / rho
            pl(1, ih) = ONE / rhosq
            pl(k - 1, j) = TWO / rho  ! Should be moved out of the loop
            pl(k - 1, ih) = -TWO / rhosq  ! Should be moved out of the loop
        else  ! XPT(J, NF) = -RHO
            d(j) = (fval(kk) + fval(k) - TWO * fbase) / rhosq
            pq(j) = (fval(kk) - fval(k)) / (TWO * rho)
            pl(1, ih) = -TWO / rhosq
            pl(k - 1, j) = HALF / rho  ! Should be moved out of the loop
            pl(k - 1, ih) = ONE / rhosq  ! Should be moved out of the loop
        end if
        pq(ih) = d(j)
        pl(k, j) = -HALF / rho
        pl(k, ih) = ONE / rhosq
    end if
end do

ih = n + 1
ip = 0
iq = 2

! Form the off-diagonal second derivatives of the initial quadratic model.
if (info == INFO_DFT) then
    do k = 2_IK * n + 2_IK, npt
        ! Pick the shift from XBASE to the next initial interpolation point that provides
        ! off-diagonal second derivatives.
        ip = ip + 1
        if (ip == iq) then
            ih = ih + 1
            ip = 1
            iq = iq + 1
        end if
        xpt(ip, k) = xw(ip)
        ! N.B.: XPT(2, NF+1) is accessed by XPT(IQ, NF+1) even if N = 1.
        xpt(iq, k) = xw(iq)
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

        ih = ih + 1
        temp = ONE / (xw(ip) * xw(iq))
        tempa = fval(k) - fbase - xw(ip) * pq(ip) - xw(iq) * pq(iq)
        ! N.B.: D(2) is accessed by D(IQ) even if N = 1.
        pq(ih) = (tempa - HALF * rhosq * (d(ip) + d(iq))) * temp
        pl(1, ih) = temp
        iw = ip + ip
        if (xw(ip) < ZERO) iw = iw + 1
        pl(iw, ih) = -temp
        iw = iq + iq
        if (xw(iq) < ZERO) iw = iw + 1
        pl(iw, ih) = -temp
        pl(k, ih) = temp
    end do
end if

nf = int(count(evaluated), kind(nf))  !!MATLAB: nf = sum(evaluated);
kopt = int(minloc(fval, mask=evaluated, dim=1), kind(kopt))
fopt = fval(kopt)
!!MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first')

end subroutine

end module initialize_mod
