module univar_mod
!--------------------------------------------------------------------------------------------------!
! This module implements functions that approximately optimizes univariate functions. They are used
! in NEWUOA (TRSAPP, BIGLAG, and BIGDEN) and BOBYQA (TRSBOX).
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
!
! Started: January 2021
!
! Last Modified: Sunday, April 24, 2022 PM02:15:25
!
! N.B.:
! 0. Here, the optimization is performed using only function values by sampling the objective
! function on a grid. Indeed, in NEWUOA and BOBYQA, the objective function optimized here is either
! a trigonometric function or a rational function. So derivatives can be used, but a simple grid
! search suffices because highly precise solutions are not necessary.
! 1. Both CIRCLE_MIN and CIRCLE_MAXABS require an input GRID_SIZE, the size of the grid used in
! the search. Powell chose GRID_SIZE = 50 in NEWUOA. MAGICALLY, this number works the best for
! NEWUOA in tests on CUTest problems. Larger (e.g., 60, 100) or smaller (e.g., 20, 40) values will
! worsen the performance of NEWUOA. Why?
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: circle_min, circle_maxabs, interval_max

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
! function at an evenly distributed "grid" on [0, 2*PI], the number of grid points being GRID_SIZE.
! Then it takes the grid point with the least value of FUN, and improves the point by a step that
! minimizes the quadratic that interpolates FUN on this point and its two nearest neighbours.
! The objective function FUN can represent the parametrization of a function defined on the circle,
! which explains the name of this function.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, PI, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : linspace
implicit none

! Inputs
procedure(FUNC_WITH_ARGS) :: fun
real(RP), intent(in) :: args(:)
integer(IK), intent(in) :: grid_size

! Outputs
real(RP) :: angle

! Local variables
character(len=*), parameter :: srname = 'CIRCLE_MIN'
integer(IK) :: k
integer(IK) :: kopt
real(RP) :: agrid(grid_size + 1)
real(RP) :: fprev
real(RP) :: fnext
real(RP) :: fopt
real(RP) :: fgrid(grid_size)
real(RP) :: step
real(RP) :: unit_angle

! Preconditions
if (DEBUGGING) then
    call assert(grid_size >= 3, 'GRID_SIZE >= 3', srname)
end if

!====================!
! Calculation starts !
!====================!

agrid = linspace(ZERO, TWO * PI, grid_size + 1_IK) ! Size: GRID_SIZE+1; the last entry will be unused
fgrid = [(fun(agrid(k), args), k=1, grid_size)]
!!MATLAB: fgrid = arrayfun(@(angle) fun(angle, args), agrid(1:grid_size));  % Same shape as `agrid`

if (all(is_nan(fgrid))) then
    angle = ZERO
    return
end if

kopt = int(minloc(fgrid, mask=(.not. is_nan(fgrid)), dim=1), IK)
fopt = fgrid(kopt)
!!MATLAB: [fopt, kopt] = min(fgrid, [], 'omitnan');
fprev = fgrid(modulo(kopt - 2_IK, grid_size) + 1)  ! Corresponds to KOPT - 1
fnext = fgrid(modulo(kopt, grid_size) + 1)  ! Corresponds to KOPT + 1

step = ZERO
if (abs(fprev - fnext) > 0) then
    fprev = fprev - fopt
    fnext = fnext - fopt
    step = HALF * (fprev - fnext) / (fprev + fnext)
end if

if (is_finite(step) .and. abs(step) > 0) then
    unit_angle = (TWO * PI) / real(grid_size, RP)
    angle = (real(kopt - 1, RP) + step) * unit_angle
    ! 1. AGRID(KOPT) = (KOPT-1) * UNIT_ANGLE. 2. ANGLE may not be in [0, 2*PI].
else
    angle = agrid(kopt)
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(is_finite(angle), 'ANGLE is finite', srname)
end if

end function circle_min


function circle_maxabs(fun, args, grid_size) result(angle)
!--------------------------------------------------------------------------------------------------!
! This function seeks an approximate maximizer of the absolute value of a 2*PI-periodic function
! FUN(X, ARGS), where the scalar X is the decision variable, and ARGS is a given vector of parameters.
! It evaluates the function at an evenly distributed "grid" on [0, 2*PI], the number of grid points
! being GRID_SIZE. Then it takes the grid point with the largest value of |FUN|, and improves the
! point by a step that maximizes the absolute value of the quadratic that interpolates FUN on this
! point and its two nearest neighbours. The objective function FUN can represent the parametrization
! of a function defined on the circle, which explains the name of this function.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, PI, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: linalg_mod, only : linspace
implicit none

! Inputs
procedure(FUNC_WITH_ARGS) :: fun
real(RP), intent(in) :: args(:)
integer(IK), intent(in) :: grid_size

! Outputs
real(RP) :: angle

! Local variables
character(len=*), parameter :: srname = 'CIRCLE_MAXABS'
integer(IK) :: k
integer(IK) :: kopt
real(RP) :: agrid(grid_size + 1)
real(RP) :: fprev
real(RP) :: fnext
real(RP) :: fopt
real(RP) :: fgrid(grid_size)
real(RP) :: step
real(RP) :: unit_angle

! Preconditions
if (DEBUGGING) then
    call assert(grid_size >= 3, 'GRID_SIZE >= 3', srname)
end if

!====================!
! Calculation starts !
!====================!

agrid = linspace(ZERO, TWO * PI, grid_size + 1_IK) ! Size: GRID_SIZE+1; the last entry is not used
fgrid = [(fun(agrid(k), args), k=1, grid_size)]
!!MATLAB: fgrid = arrayfun(@(angle) fun(angle, args), agrid(1:grid_size));  % Same shape as `agrid`

if (all(is_nan(fgrid))) then
    angle = ZERO
    return
end if

kopt = int(maxloc(abs(fgrid), mask=(.not. is_nan(fgrid)), dim=1), IK)
!!MATLAB: [~, kopt] = max(abs(fgrid), [], 'omitnan');
fopt = fgrid(kopt)
fprev = fgrid(modulo(kopt - 2_IK, grid_size) + 1)  ! Corresponds to KOPT - 1
fnext = fgrid(modulo(kopt, grid_size) + 1)  ! Corresponds to KOPT + 1

step = ZERO
if (abs(fprev - fnext) > 0) then
    fprev = fprev - fopt
    fnext = fnext - fopt
    step = HALF * (fprev - fnext) / (fprev + fnext)
end if

if (is_finite(step) .and. abs(step) > 0) then
    unit_angle = (TWO * PI) / real(grid_size, RP)
    angle = (real(kopt - 1, RP) + step) * unit_angle
    ! 1. AGRID(KOPT) = (KOPT-1) * UNIT_ANGLE. 2. ANGLE may not be in [0, 2*PI].
else
    angle = agrid(kopt)
end if

! Postconditions
if (DEBUGGING) then
    call assert(is_finite(angle), 'ANGLE is finite', srname)
end if

end function circle_maxabs


function interval_max(fun, lb, ub, args, grid_size) result(x)
!--------------------------------------------------------------------------------------------------!
! This function seeks an approximate maximizer of a function F(X, ARGS) for X in [LB, UB], where
! ARGS is a vector of parameters. It evaluates the function at an evenly distributed "grid" on
! [LB, UB], the number of grid points being GRID_SIZE. Then it takes the grid point with the largest
! value of FUN, and improves the point by a step that maximizes the quadratic that interpolates
! FUN on this point and its two nearest neighbours unless the point is LB or UB.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, HALF, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
use, non_intrinsic :: linalg_mod, only : linspace
implicit none

! Inputs
procedure(FUNC_WITH_ARGS) :: fun
real(RP), intent(in) :: lb
real(RP), intent(in) :: ub
real(RP), intent(in) :: args(:)
integer(IK), intent(in) :: grid_size

! Outputs
real(RP) :: x

! Local variables
character(len=*), parameter :: srname = 'INTERVAL_MAX'
integer(IK) :: k
integer(IK) :: kopt
real(RP) :: fgrid(grid_size)
real(RP) :: fopt
real(RP) :: fnext
real(RP) :: fprev
real(RP) :: step
real(RP) :: xgrid(grid_size)

! Preconditions
if (DEBUGGING) then
    call assert(is_finite(lb) .and. is_finite(ub) .and. lb <= ub, 'LB <= UB and they are finite', srname)
    call assert(grid_size >= 3, 'GRID_SIZE >= 3', srname)
end if

!====================!
! Calculation starts !
!====================!

if (ub <= lb) then
    x = lb
    return
end if

xgrid = linspace(lb, ub, grid_size)
fgrid = [(fun(xgrid(k), args), k=1, grid_size)]
!!MATLAB: fgrid = arrayfun(@(x) fun(x, args), xgrid(1:grid_size));  % Same shape as `xgrid`

if (all(is_nan(fgrid))) then
    x = lb
    return
end if

kopt = int(maxloc(fgrid, mask=(.not. is_nan(fgrid)), dim=1), IK)
fopt = fgrid(kopt)
!!MATLAB: [fopt, kopt] = min(fgrid, [], 'omitnan');

if (kopt == 1) then
    x = lb
elseif (kopt == grid_size) then
    x = ub
else
    fprev = fgrid(kopt - 1)
    fnext = fgrid(kopt + 1)
    step = ZERO
    if (abs(fprev - fnext) > 0) then
        step = HALF * ((fnext - fprev) / (fopt + fopt - fprev - fnext))
    end if
    if (is_finite(step) .and. abs(step) > 0) then
        x = lb + (ub - lb) * (real(kopt - 1, RP) + step) / real(grid_size - 1, RP)
        ! N.B.: 1. XGRID(KOPT) = LB + (UB-LB)*(KOPT - 1)/(GRID_SIZE -1)
        ! 2. XGRID(KOPT-1) <= X <= XGRID(KOPT+1), as X maximizes the quadratic interpolant.
    else
        x = xgrid(kopt)
    end if
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(lb <= x .and. x <= ub, 'LB <= X <= UB', srname)
end if

end function interval_max


end module univar_mod
