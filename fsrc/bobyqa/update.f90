module update_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the update of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Saturday, March 05, 2022 PM04:51:13
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: update


contains


subroutine update(knew, beta, denom, bmat, vlag, zmat)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert

implicit none

! Inputs
integer(IK), intent(in) :: knew
real(RP), intent(in) :: beta
real(RP), intent(in) :: denom

! In-outputs
real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(inout) :: vlag(:)  ! VLAG(NPT + N)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Local variables
character(len=*), parameter :: srname = 'UPDATE'
integer(IK) :: n
integer(IK) :: npt
real(RP) :: w(size(vlag))
real(RP) :: alpha, tau, temp, tempa, tempb, ztest
integer(IK) :: i, j, jp, k

! Sizes.
n = int(size(bmat, 1), kind(n))
npt = int(size(bmat, 2) - size(bmat, 1), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(vlag) == npt + n, 'SIZE(VLAG) == NPT + N', srname)
end if



!     The arrays BMAT and ZMAT are updated, as required by the new position
!     of the interpolation point that has the index KNEW. The vector VLAG has
!     N+NPT components, set on entry to the first NPT and last N components
!     of the product Hw in equation (4.11) of the Powell (2006) paper on
!     NEWUOA. Further, BETA is set on entry to the value of the parameter
!     with that name, and DENOM is set to the denominator of the updating
!     formula. Elements of ZMAT may be treated as zero if their moduli are
!     at most ZTEST. The first NDIM elements of W are used for working space.
!
!     Set some constants.
!
ztest = ZERO
do k = 1, npt
    do j = 1, npt - n - 1_IK
        ztest = max(ztest, abs(zmat(k, j)))
    end do
end do
ztest = 1.0E-20_RP * ztest
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
! Zaikun 2019-08-15: JL is never used
!      JL=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j = 2, npt - n - 1_IK
    if (abs(zmat(knew, j)) > ztest) then
        temp = sqrt(zmat(knew, 1)**2 + zmat(knew, j)**2)
        tempa = zmat(knew, 1) / temp
        tempb = zmat(knew, j) / temp
        do i = 1, npt
            temp = tempa * zmat(i, 1) + tempb * zmat(i, j)
            zmat(i, j) = tempa * zmat(i, j) - tempb * zmat(i, 1)
            zmat(i, 1) = temp
        end do
    end if
    zmat(knew, j) = ZERO
end do
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
do i = 1, npt
    w(i) = zmat(knew, 1) * zmat(i, 1)
end do
alpha = w(knew)
tau = vlag(knew)
vlag(knew) = vlag(knew) - ONE
!
!     Complete the updating of ZMAT.
!
temp = sqrt(denom)
tempb = zmat(knew, 1) / temp
tempa = tau / temp
do i = 1, npt
    zmat(i, 1) = tempa * zmat(i, 1) - tempb * vlag(i)
end do
!
!     Finally, update the matrix BMAT.
!
do j = 1, n
    jp = npt + j
    w(jp) = bmat(j, knew)
    tempa = (alpha * vlag(jp) - tau * w(jp)) / denom
    tempb = (-beta * w(jp) - tau * vlag(jp)) / denom
    do i = 1, jp
        bmat(j, i) = bmat(j, i) + tempa * vlag(i) + tempb * w(i)
        if (i > npt) bmat(i - npt, jp) = bmat(j, i)
    end do
end do
return

end subroutine update


end module update_mod
