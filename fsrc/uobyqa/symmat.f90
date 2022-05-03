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
! Last Modified: Wednesday, May 04, 2022 AM12:43:47
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: vec2mat, mat2vec


contains


function vec2mat(vec) result(mat)
!--------------------------------------------------------------------------------------------------!
! This module transforms a vector to a symmetric matrix with the vector storing the upper triangular
! part of the matrix column by column.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : issymmetric
implicit none
! Inputs
real(RP), intent(in) :: vec(:)
! Outputs
real(RP) :: mat((floor(sqrt(real(8 * size(vec) + 1))) - 1) / 2, (floor(sqrt(real(8 * size(vec) + 1))) - 1) / 2)
! Local variables
character(len=*), parameter :: srname = 'MAT2VEC'
integer(IK) :: i
integer(IK) :: ih
integer(IK) :: j
integer(IK) :: n

! Sizes
n = int(size(mat, 1), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(size(vec) == n * (n + 1) / 2, 'SIZE(VEC) = N*(N+1)/2', 'vec2mat')
end if

!====================!
! Calculation starts !
!====================!

ih = 0_IK
do j = 1_IK, n
    do i = 1_IK, j
        ih = ih + 1_IK
        mat(i, j) = vec(ih)
        mat(j, i) = vec(ih)
    end do
end do

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(issymmetric(mat), 'MAT is symmetric', srname)
end if
end function vec2mat


function mat2vec(mat) result(vec)
!--------------------------------------------------------------------------------------------------!
! This module transforms a symmetric matrix to a vector that stores the upper triangular part of the
! matrix column by column.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : issymmetric
implicit none
! Inputs
real(RP), intent(in) :: mat(:, :)
! Outputs
real(RP) :: vec((size(mat, 1) * (size(mat, 1) + 1)) / 2)
! Local variables
character(len=*), parameter :: srname = 'VEC2MAT'
integer(IK) :: i
integer(IK) :: ih
integer(IK) :: n
integer(IK) :: j

! Preconditions
if (DEBUGGING) then
    call assert(issymmetric(mat), 'MAT is symmetric', srname)
end if

!====================!
! Calculation starts !
!====================!

n = int(size(mat, 1), kind(n))
ih = 0_IK
do j = 1_IK, n
    do i = 1_IK, j
        ih = ih + 1_IK
        vec(ih) = mat(i, j)
    end do
end do

!====================!
! Calculation starts !
!====================!

end function mat2vec

end module symmat_mod
