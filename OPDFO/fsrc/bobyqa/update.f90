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
! Last Modified: Wednesday, April 13, 2022 AM01:59:12
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: update


contains


subroutine update(n, npt, bmat, zmat, ndim, vlag, beta, denom, knew, w)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO
use, non_intrinsic :: linalg_mod, only : planerot

implicit none

integer(IK), intent(in) :: knew
integer(IK), intent(in) :: n
integer(IK), intent(in) :: ndim
integer(IK), intent(in) :: npt
real(RP), intent(in) :: beta
real(RP), intent(in) :: denom
real(RP), intent(inout) :: bmat(n, npt + n)
real(RP), intent(inout) :: vlag(npt + n)
real(RP), intent(inout) :: w(npt + n)
real(RP), intent(inout) :: zmat(npt, npt - n - 1_IK)

! local variables
real(RP) :: alpha, tau, temp, tempa, tempb, ztest, grot(2, 2), zmatk1
integer(IK) :: i, j, jp, k, nptm


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
nptm = npt - n - 1
ztest = ZERO
do k = 1, npt
    do j = 1, nptm
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
do j = 2, nptm
    if (abs(zmat(knew, j)) > ztest) then
        !temp = sqrt(zmat(knew, 1)**2 + zmat(knew, j)**2)
        !tempa = zmat(knew, 1) / temp
        !tempb = zmat(knew, j) / temp
        grot = planerot(zmat(knew, [1, j]))
        tempa = grot(1, 1)
        tempb = grot(1, 2)
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
zmatk1 = zmat(knew, 1)
do i = 1, npt
    !zmat(i, 1) = tempa * zmat(i, 1) - tempb * vlag(i)
    zmat(i, 1) = (tau * zmat(i, 1) - zmatk1 * vlag(i)) / temp
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
