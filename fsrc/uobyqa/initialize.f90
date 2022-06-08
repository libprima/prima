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
! Last Modified: Wednesday, June 08, 2022 AM08:44:19
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
real(RP) :: f
real(RP) :: x(size(x0))

integer(IK) :: iw, ih, ip, iq, i, j
real(RP) :: rho, rhosq, fplus, fbase, xw(size(x0)), d(size(x0)), temp, tempa

n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = max(maxxhist, maxfhist)

rho = rhobeg
rhosq = rho**2

! Initialization. NF is the number of function calculations so far. The least function value so far,
! the corresponding X and its index are noted in FOPT, XOPT, and KOPT respectively.
info = INFO_DFT
xbase = x0
xpt = ZERO
pl = ZERO
j = 0
ih = n
! In the following loop, FPLUS is set to F(X + RHO*e_I) when NF = 2*I, and the value of FPLUS is
! used subsequently when NF = 2*I + 1.
fplus = ZERO ! This initial value is not used but to entertain the Fortran compilers.
do nf = 1, 2_IK * n + 1_IK
    ! Pick the shift from XBASE to the next initial interpolation point that provides diagonal
    ! second derivatives.
    if (nf > 1) then
        if (modulo(nf, 2_IK) == 1_IK) then
            if (fplus < fbase) then
                xw(j) = rho
                xpt(j, nf) = TWO * rho
            else
                xw(j) = -rho
                xpt(j, nf) = -rho
            end if
        elseif (j < n) then
            j = j + 1
            xpt(j, nf) = rho
        end if
    end if
    x = xbase + xpt(:, nf)

    if (is_nan(sum(abs(x)))) then
        f = sum(x) ! Set F to NaN
        if (nf == 1) then
            kopt = 1
            fopt = f
        end if
        info = NAN_INF_X
        exit
    end if
    call evaluate(calfun, x, f)
    call savehist(nf, x, xhist, f, fhist)
    if (is_nan(f) .or. is_posinf(f)) then
        if (nf == 1) then
            kopt = 1
            fopt = f
        end if
        info = NAN_INF_F
        exit
    end if

    if (f <= ftarget) then
        info = FTARGET_ACHIEVED
        kopt = nf
        fopt = f
        exit
    end if

    if (nf == 1) then
        fopt = f
        kopt = nf
        fbase = f
    elseif (f < fopt) then
        fopt = f
        kopt = nf
    end if

    if (nf >= maxfun) then
        info = MAXFUN_REACHED
        exit
    end if

    if (modulo(nf, 2_IK) == 0_IK) then
        fplus = f  ! FPLUS = F(X + RHO*e_I) with I = NF/2.
        cycle
    end if

    ! Form the gradient and diagonal second derivatives of the quadratic model and Lagrange functions.
    if (j >= 1 .and. nf >= 3) then  ! NF >= 3 is implied by J >= 1. We prefer to impose it explicitly.
        ih = ih + j
        if (xpt(j, nf) > 0) then  ! XPT(J, NF) = 2*RHO
            pq(j) = (4.0_RP * fplus - 3.0_RP * fbase - f) / (TWO * rho)
            d(j) = (fbase + f - TWO * fplus) / rhosq
            pl(1, j) = -1.5_RP / rho
            pl(1, ih) = ONE / rhosq
            pl(nf - 1, j) = TWO / rho  ! Should be moved out of the loop
            pl(nf - 1, ih) = -TWO / rhosq  ! Should be moved out of the loop
        else  ! XPT(J, NF) = -RHO
            d(j) = (fplus + f - TWO * fbase) / rhosq
            pq(j) = (fplus - f) / (TWO * rho)
            pl(1, ih) = -TWO / rhosq
            pl(nf - 1, j) = HALF / rho  ! Should be moved out of the loop
            pl(nf - 1, ih) = ONE / rhosq  ! Should be moved out of the loop
        end if
        pq(ih) = d(j)
        pl(nf, j) = -HALF / rho
        pl(nf, ih) = ONE / rhosq
    end if
end do

ih = n + 1
ip = 0
iq = 2

! Form the off-diagonal second derivatives of the initial quadratic model.
if (info == INFO_DFT) then
    do nf = 2_IK * n + 2_IK, npt
        ! Pick the shift from XBASE to the next initial interpolation point that provides
        ! off-diagonal second derivatives.
        ip = ip + 1
        if (ip == iq) then
            ih = ih + 1
            ip = 1
            iq = iq + 1
        end if
        xpt(ip, nf) = xw(ip)
        ! N.B.: XPT(2, NF+1) is accessed by XPT(IQ, NF+1) even if N = 1.
        xpt(iq, nf) = xw(iq)
        x = xbase + xpt(:, nf)

        if (is_nan(sum(abs(x)))) then
            f = sum(x) ! Set F to NaN
            info = NAN_INF_X
            exit
        end if

        call evaluate(calfun, x, f)
        call savehist(nf, x, xhist, f, fhist)
        if (is_nan(f) .or. is_posinf(f)) then
            info = NAN_INF_F
            exit
        end if

        if (f <= ftarget) then
            info = FTARGET_ACHIEVED
            kopt = nf
            fopt = f
            exit
        end if

        if (f < fopt) then
            fopt = f
            kopt = nf
        end if
        if (nf >= maxfun) then
            info = MAXFUN_REACHED
            exit
        end if

        ih = ih + 1
        temp = ONE / (xw(ip) * xw(iq))
        tempa = f - fbase - xw(ip) * pq(ip) - xw(iq) * pq(iq)
        ! N.B.: D(2) is accessed by D(IQ) even if N = 1.
        pq(ih) = (tempa - HALF * rhosq * (d(ip) + d(iq))) * temp
        pl(1, ih) = temp
        iw = ip + ip
        if (xw(ip) < ZERO) iw = iw + 1
        pl(iw, ih) = -temp
        iw = iq + iq
        if (xw(iq) < ZERO) iw = iw + 1
        pl(iw, ih) = -temp
        pl(nf, ih) = temp
    end do
end if

!--------------------------------------------------------------------------------------------------!
! When the loop exits, the value of NF is not specified by the standard. With gfortran, it will be
! NPT+1, which is not proper for the subsequent use. !!!
nf = min(nf, npt) !!!

end subroutine

end module initialize_mod
