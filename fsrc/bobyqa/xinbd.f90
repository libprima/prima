module xinbd_mod

implicit none

private
public :: xinbd


contains


function xinbd(xbase, step, xl, xu, sl, su) result(x)
!--------------------------------------------------------------------------------------------------!
! This function sets X to XBASE + STEP, paying careful attention to the following bounds.
! 1. XBASE is a point between XL and XU (guaranteed);
! 2. STEP is a step between SL and SU (may be with rounding errors);
! 3. SL = XL - XBASE, SU = XU - XBASE;
! 4. X should be between XL and XU.
!--------------------------------------------------------------------------------------------------!
! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, EPS, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : trueloc
use, non_intrinsic :: infnan_mod, only : is_finite

implicit none

! Inputs
real(RP), intent(in) :: xbase(:)
real(RP), intent(in) :: step(:)
real(RP), intent(in) :: xl(:)
real(RP), intent(in) :: xu(:)
real(RP), intent(in) :: sl(:)
real(RP), intent(in) :: su(:)

! Outputs
real(RP) :: x(size(xbase))

! Local variables
character(len=*), parameter :: srname = 'XINBD'
integer(IK) :: n
real(RP) :: s(size(xbase))

! Sizes
n = int(size(xbase), kind(n))

! Preconditions
if (DEBUGGING) then
    call assert(all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == N == SIZE(XU)', srname)
    call assert(all(xbase >= xl .and. xbase <= xu), 'XL <= XBASE <= XU', srname)
    call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(all(step + 1.0E2_RP * EPS * max(ONE, abs(step)) >= sl .and. &
        & step - 1.0E2_RP * EPS * max(ONE, abs(step)) <= su), 'SL <= STEP <= SU', srname)
end if

!====================!
! Calculation starts !
!====================!

s = max(sl, min(su, step))
x = max(xl, min(xu, xbase + s))
x(trueloc(s <= sl)) = xl(trueloc(s <= sl))
x(trueloc(s >= su)) = xu(trueloc(s >= su))

!====================!
!  Calculation ends  !
!====================!

if (DEBUGGING) then
    call assert(size(x) == n .and. all(x >= xl .and. x <= xu), 'SIZE(X) == N, XL <= X <= XU', srname)
    call assert(all(x <= xl .or. step > sl), 'X == XL if STEP <= SL', srname)
    call assert(all(x >= xu .or. step < su), 'X == XU if STEP >= SU', srname)
end if

end function xinbd


end module xinbd_mod
