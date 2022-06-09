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
! Last Modified: Thursday, June 09, 2022 PM12:15:03
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: initxf, initq, initl


contains


subroutine initxf(calfun, iprint, maxfun, ftarget, rhobeg, x0, kopt, nf, fhist, fval, xbase, &
    & xhist, xpt, info)

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: info_mod, only : INFO_DFT
use, non_intrinsic :: linalg_mod, only : eye, trueloc, linspace
use, non_intrinsic :: output_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ

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
real(RP), intent(out) :: fval(:)
real(RP), intent(out) :: xbase(:)   ! XBASE(N)
real(RP), intent(out) :: xhist(:, :)    ! XHIST(N, MAXXHIST)
real(RP), intent(out) :: xpt(:, :)  ! XPT(N, NPT)

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

integer(IK) :: ip, iq, kk(size(x0))
real(RP) :: rho, xw(size(x0))

n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = max(maxxhist, maxfhist)

rho = rhobeg

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

! Set XPT(:, 1 : 2*N+1) and FVAL(:, 1 : 2*N+1).
xpt = ZERO
kk = linspace(2_IK, 2_IK * n, n); 
xpt(:, kk) = rho * eye(n)
do k = 1, 2_IK * n + 1_IK
    if (k >= 3 .and. modulo(k, 2_IK) == 1) then
        if (fval(k - 1) < fval(1)) then
            xpt((k - 1) / 2, k) = TWO * rho
        else
            xpt((k - 1) / 2, k) = -rho
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
end do

if (info == INFO_DFT .and. n > 1) then
    xw = -rho
    xw(trueloc(fval(kk) < fval(1))) = rho

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
!!MATLAB: fopt = min(fval(evaluated)); kopt = find(evaluated & ~(fval > fopt), 1, 'first')

end subroutine initxf


subroutine initq(fval, xpt, pq, info)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_MODEL

implicit none

! Inputs
real(RP), intent(in) :: fval(:)  ! XPT(N, NPT)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! Outputs
integer(IK), intent(out), optional :: info
real(RP), intent(out) :: pq(:)  ! PQ((N + 1) * (N + 2) / 2 - 1)

! Local variables
character(len=*), parameter :: srname = 'INITL'
integer(IK) :: n
integer(IK) :: npt

integer(IK) :: ih, ip, iq, k0, k1, k
real(RP) :: rho, rhosq, fbase, deriv(size(xpt, 1))

n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

rho = maxval(abs(xpt(:, 2)))
rhosq = rho**2
fbase = fval(1)

! Form the gradient and diagonal second derivatives of the quadratic model.
do k = 1, n
    k0 = 2_IK * k
    k1 = 2_IK * k + 1_IK
    ih = n + k * (k + 1_IK) / 2_IK
    if (xpt(k, k1) > 0) then  ! XPT(K, K1) = 2*RHO
        deriv(k) = (fbase + fval(k1) - TWO * fval(k0)) / rhosq
        pq(k) = (4.0_RP * fval(k0) - 3.0_RP * fbase - fval(k1)) / (TWO * rho)
    else  ! XPT(K, K1) = -RHO
        deriv(k) = (fval(k0) + fval(k1) - TWO * fbase) / rhosq
        pq(k) = (fval(k0) - fval(k1)) / (TWO * rho)
    end if
    pq(ih) = deriv(k)
end do

! Form the off-diagonal second derivatives of the initial quadratic model.
if (n > 1) then
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
        pq(ih) = (fval(k) - fbase - xpt(ip, k) * pq(ip) - xpt(iq, k) * pq(iq) &
            & - HALF * rhosq * (deriv(ip) + deriv(iq))) / (xpt(ip, k) * xpt(iq, k))
    end do
end if

if (present(info)) then
    if (is_nan(sum(abs(pq)))) then
        info = NAN_MODEL
    else
        info = INFO_DFT
    end if
end if

end subroutine initq


subroutine initl(xpt, pl, info)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
use, non_intrinsic :: info_mod, only : INFO_DFT, NAN_MODEL

implicit none

! Inputs
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! Outputs
integer(IK), intent(out), optional :: info
real(RP), intent(out) :: pl(:, :)  ! PL((N + 1) * (N + 2) / 2, (N + 1) * (N + 2) / 2 - 1)

! Local variables
character(len=*), parameter :: srname = 'INITL'
integer(IK) :: k
integer(IK) :: n
integer(IK) :: npt

integer(IK) :: ih, ip, iq, k0, k1
real(RP) :: rho, rhosq, temp

n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))

rho = maxval(abs(xpt(:, 2)))
rhosq = rho**2

pl = ZERO

! Form the gradient and diagonal second derivatives of the Lagrange functions.
do k = 1, n
    k0 = 2_IK * k
    k1 = 2_IK * k + 1_IK
    ih = n + k * (k + 1_IK) / 2_IK
    if (xpt(k, k1) > 0) then  ! XPT(K, K1) = 2*RHO
        pl(1, k) = -1.5_RP / rho
        pl(1, ih) = ONE / rhosq
        pl(k0, k) = TWO / rho
        pl(k0, ih) = -TWO / rhosq
    else  ! XPT(K, K1) = -RHO
        pl(1, ih) = -TWO / rhosq
        pl(k0, k) = HALF / rho
        pl(k0, ih) = ONE / rhosq
    end if
    pl(k1, k) = -HALF / rho
    pl(k1, ih) = ONE / rhosq
end do

! Form the off-diagonal second derivatives of the Lagrange functions.
if (n > 1) then
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
        temp = ONE / (xpt(ip, k) * xpt(iq, k))
        pl(1, ih) = temp
        pl(k, ih) = temp

        if (xpt(ip, k) < 0) then
            pl(2_IK * ip + 1_IK, ih) = -temp
        else
            pl(2_IK * ip, ih) = -temp
        end if

        if (xpt(iq, k) < 0) then
            pl(2_IK * iq + 1_IK, ih) = -temp
        else
            pl(2_IK * iq, ih) = -temp
        end if
    end do
end if

if (present(info)) then
    if (is_nan(sum(abs(pl)))) then
        info = NAN_MODEL
    else
        info = INFO_DFT
    end if
end if

end subroutine initl

end module initialize_mod
