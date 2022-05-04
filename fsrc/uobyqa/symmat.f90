module symmat_mod
!--------------------------------------------------------------------------------------------------!
! This module provides functions that transforms between a vector a symmetric matrix with the vector
! storing the upper triangular part of the matrix column by column.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the UOBYQA paper.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Wednesday, May 04, 2022 PM04:46:00
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: vec2smat, smat2vec, smat_mul_vec


contains


function vec2smat(vec) result(smat)
!--------------------------------------------------------------------------------------------------!
! This function transforms a vector VEC to a symmetric matrix SMAT with the vector storing the upper
! triangular part of the matrix column by column.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : issymmetric
implicit none
! Inputs
real(RP), intent(in) :: vec(:)
! Outputs
real(RP) :: smat((floor(sqrt(real(8 * size(vec) + 1))) - 1) / 2, (floor(sqrt(real(8 * size(vec) + 1))) - 1) / 2)
! Local variables
character(len=*), parameter :: srname = 'SMAT2VEC'
integer(IK) :: ih
integer(IK) :: j
integer(IK) :: n

! Sizes
n = int(size(smat, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(vec) == n * (n + 1) / 2, 'SIZE(VEC) = N*(N+1)/2', srname)
end if

!====================!
! Calculation starts !
!====================!

do j = 1_IK, n
    ih = (j - 1_IK) * j / 2_IK
    smat(1:j, j) = vec(ih + 1:ih + j)
    smat(j, 1:j - 1) = smat(1:j - 1, j)
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(issymmetric(smat), 'SMAT is symmetric', srname)
end if
end function vec2smat


function smat2vec(smat) result(vec)
!--------------------------------------------------------------------------------------------------!
! This function transforms a symmetric matrix SMAT to a vector VEC that stores the upper triangular
! part of the matrix column by column.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : issymmetric
implicit none
! Inputs
real(RP), intent(in) :: smat(:, :)
! Outputs
real(RP) :: vec((size(smat, 1) * (size(smat, 1) + 1)) / 2)
! Local variables
character(len=*), parameter :: srname = 'SMAT2VEC'
integer(IK) :: ih
integer(IK) :: n
integer(IK) :: j

! Preconditions
if (DEBUGGING) then
    call assert(issymmetric(smat), 'SMAT is symmetric', srname)
end if

!====================!
! Calculation starts !
!====================!

n = int(size(smat, 1), kind(n))
do j = 1_IK, n
    ih = (j - 1_IK) * j / 2_IK
    vec(ih + 1:ih + j) = smat(1:j, j)
end do

!====================!
! Calculation starts !
!====================!

end function smat2vec


function smat_mul_vec(smatv, x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function calculates the product of a symmetric matrix and a vector X,  with the upper
! triangular part of the matrix stored in the vector SMATV column by column.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : inprod
implicit none
! Inputs
real(RP), intent(in) :: smatv(:)
real(RP), intent(in) :: x(:)
! Outputs
real(RP) :: y(size(x))
! Local variables
character(len=*), parameter :: srname = 'SMAT_MUL_VEC'
integer(IK) :: ih
integer(IK) :: n
integer(IK) :: j

! Sizes
n = int(size(x), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(smatv) == n * (n + 1_IK) / 2_IK, 'SIZE(SMATV) = N*(N+1)/2', srname)
end if

!====================!
! Calculation starts !
!====================!

do j = 1, n
    ih = (j - 1_IK) * j / 2_IK
    y(j) = inprod(smatv(ih + 1:ih + j), x(1:j))
    y(1:j - 1) = y(1:j - 1) + x(j) * smatv(ih + 1:ih + j - 1)
end do

!====================!
! Calculation starts !
!====================!

end function smat_mul_vec


end module symmat_mod
