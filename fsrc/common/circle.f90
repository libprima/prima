module circle_mod
!--------------------------------------------------------------------------------------------------!
! This module implements functions that approximately optimizes functions on circles (equivalently,
! 2*PI-periodic functions). They are used in TRSAPP, BIGLAG, and BIGDEN of NEWUOA.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: January 2021
!
! Last Modified: Thursday, January 06, 2022 AM12:01:38
!
! N.B.: Both CIRCLE_MIN and CIRCLE_MAXABS require an input GRID_SIZE, the size of the grid used in
! the search. Powell chose GRID_SIZE = 50 in NEWUOA. MAGICALLY, this number works the best for
! NEWUOA in tests on CUTest problems. Larger (e.g., 60, 100) or smaller (e.g., 20, 40) values will
! worsen the performance of NEWUOA. Why?
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: circle_min, circle_maxabs

abstract interface
    function FUNC_WITH_ARGS(x, args) result(f)
    use, non_intrinsic :: consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x
    real(RP), intent(in) :: args(:)
    real(RP) :: f
    end function FUNC_WITH_ARGS
end interface


contains


function circle_min(fun, args, grid_size) result(angle)
!--------------------------------------------------------------------------------------------------!
! This function seeks an approximate minimizer of a 2*PI-periodic function FUN(X, ARGS), where the
! scalar X is the decision variable, and ARGS is a given vector of parameters. It evaluates the
! function at an evenly distributed "grid" on [0, 2*PI], GRID_SIZE being the number of grid points.
! Then it takes the grid point with the least value of FUN, and improves the point by a step that
! minimizes the quadratic that interpolates FUN on this point and its two nearest neighbours.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, PI, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
procedure(FUNC_WITH_ARGS) :: fun
real(RP), intent(in) :: args(:)
integer(IK), intent(in) :: grid_size

! Outputs
real(RP) :: angle

! Local variables
character(len=*), parameter :: srname = 'CIRCLE_MIN'
integer(IK) :: i
integer(IK) :: iopt
real(RP) :: angles(grid_size)
real(RP) :: fa
real(RP) :: fb
real(RP) :: fopt
real(RP) :: fval(grid_size)
real(RP) :: step
real(RP) :: unit_angle

! Preconditions
if (DEBUGGING) then
    call assert(grid_size >= 3, 'GRID_SIZE >= 3', srname)
end if

!====================!
! Calculation starts !
!====================!

unit_angle = (TWO * PI) / real(grid_size, RP)
angles = unit_angle * real([(i, i=1, grid_size)] - 1, RP)
fval = [(fun(angles(i), args), i=1, grid_size)]

if (all(is_nan(fval))) then
    angle = ZERO
    return
end if

iopt = int(minloc(fval, mask=(.not. is_nan(fval)), dim=1) - 1, IK)
fopt = fval(iopt + 1)
fa = fval(modulo(iopt - 1_IK, grid_size) + 1)
fb = fval(modulo(iopt + 1_IK, grid_size) + 1)
if (abs(fa - fb) > ZERO) then
    fa = fa - fopt
    fb = fb - fopt
    step = HALF * (fa - fb) / (fa + fb)
else
    step = ZERO
end if
angle = unit_angle * (real(iopt, RP) + step)  ! It may not be in [0, 2*PI].

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(.not. is_nan(angle), 'ANGLE is not NaN', srname)
end if

end function circle_min


function circle_maxabs(fun, args, grid_size) result(angle)
!--------------------------------------------------------------------------------------------------!
! This function seeks an approximate maximizer of the absolute value of a 2*PI-periodic function
! FUN(X, ARGS), where the scalar X is the decision variable, and ARGS is a given vector of parameters.
! It evaluates the function at an evenly distributed "grid" on [0, 2*PI], GRID_SIZE being the number
! of grid points. It takes the grid point with the largest value of |FUN|, and then improves the
! point by a step that maximizes the absolute value of the quadratic that interpolates FUN on this
! point and its two nearest neighbours.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, PI, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
procedure(FUNC_WITH_ARGS) :: fun
real(RP), intent(in) :: args(:)
integer(IK), intent(in) :: grid_size

! Outputs
real(RP) :: angle

! Local variables
character(len=*), parameter :: srname = 'CIRCLE_MAXABS'
integer(IK) :: i
integer(IK) :: iopt
real(RP) :: angles(grid_size)
real(RP) :: fa
real(RP) :: fb
real(RP) :: fopt
real(RP) :: fval(grid_size)
real(RP) :: step
real(RP) :: unit_angle

! Preconditions
if (DEBUGGING) then
    call assert(grid_size >= 3, 'GRID_SIZE >= 3', srname)
end if

!====================!
! Calculation starts !
!====================!

unit_angle = (TWO * PI) / real(grid_size, RP)
angles = unit_angle * real([(i, i=1, grid_size)] - 1, RP)
fval = [(fun(angles(i), args), i=1, grid_size)]

if (all(is_nan(fval))) then
    angle = ZERO
    return
end if

iopt = int(maxloc(abs(fval), mask=(.not. is_nan(fval)), dim=1) - 1, IK)
fopt = fval(iopt + 1)
fa = fval(modulo(iopt - 1_IK, grid_size) + 1)
fb = fval(modulo(iopt + 1_IK, grid_size) + 1)
if (abs(fa - fb) > ZERO) then
    fa = fa - fopt
    fb = fb - fopt
    step = HALF * (fa - fb) / (fa + fb)
else
    step = ZERO
end if
angle = unit_angle * (real(iopt, RP) + step)  ! It may not be in [0, 2*PI].

! Postconditions
if (DEBUGGING) then
    call assert(.not. is_nan(angle), 'ANGLE is not NaN', srname)
end if

end function circle_maxabs


end module circle_mod
