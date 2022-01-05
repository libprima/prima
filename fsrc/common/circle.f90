module circle_mod
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
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, PI!, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
procedure(FUNC_WITH_ARGS) :: fun
real(RP), intent(in) :: args(:)
integer(IK), intent(in) :: grid_size

! Output
real(RP) :: angle

! Local variables
!character(len=*), parameter :: srname = 'CIRCLE_MIN'
integer(IK) :: i
integer(IK) :: iopt
real(RP) :: angles(grid_size)
real(RP) :: fa
real(RP) :: fb
real(RP) :: fopt
real(RP) :: fval(grid_size)
real(RP) :: step
real(RP) :: unit_angle

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
end function circle_min


function circle_maxabs(fun, args, grid_size) result(angle)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, PI!, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
procedure(FUNC_WITH_ARGS) :: fun
real(RP), intent(in) :: args(:)
integer(IK), intent(in) :: grid_size

! Output
real(RP) :: angle

! Local variables
!character(len=*), parameter :: srname = 'CIRCLE_MAXABS'
integer(IK) :: i
integer(IK) :: iopt
real(RP) :: angles(grid_size)
real(RP) :: fa
real(RP) :: fb
real(RP) :: fopt
real(RP) :: fval(grid_size)
real(RP) :: step
real(RP) :: unit_angle

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
end function circle_maxabs

end module circle_mod
