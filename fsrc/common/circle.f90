module circle_mod
implicit none
private
public :: circle_search

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


function circle_search(fun, args, gridsize) result(angle)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, TWO, HALF, PI!, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan
implicit none

! Inputs
procedure(FUNC_WITH_ARGS) :: fun
real(RP), intent(in) :: args(:)
integer(IK), intent(in) :: gridsize

! Output
real(RP) :: angle

! Local variables
!character(len=*), parameter :: srname = 'CIRCLE_SEARCH'
integer(IK) :: i
integer(IK) :: imin
real(RP) :: angles(gridsize)
real(RP) :: funa
real(RP) :: funb
real(RP) :: funmin
real(RP) :: funvals(gridsize)
real(RP) :: istep
real(RP) :: unitang

unitang = (TWO * PI) / real(gridsize, RP)
angles = unitang * real([(i, i=1, gridsize)] - 1, RP)
funvals = [(fun(angles(i), args), i=1, gridsize)]
if (all(is_nan(funvals))) then
    angle = ZERO
    return
end if
imin = int(minloc(funvals, mask=(.not. is_nan(funvals)), dim=1) - 1, IK)
funmin = funvals(imin + 1)
funa = funvals(modulo(imin - 1_IK, gridsize) + 1)
funb = funvals(modulo(imin + 1_IK, gridsize) + 1)
if (abs(funa - funb) > ZERO) then
    funa = funa - funmin
    funb = funb - funmin
    istep = HALF * (funa - funb) / (funa + funb)
else
    istep = ZERO
end if
angle = unitang * (real(imin, RP) + istep)
end function circle_search


end module circle_mod
