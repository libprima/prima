subroutine construct_fletcheq2(prob)
use, non_intrinsic :: consts_mod, only : ONE, HALF, HUGENUM
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
integer(IK) :: n

prob % probname = 'fletcheq2'
prob % probtype = 'n'
prob % m = 3
prob % n = 3
call safealloc(prob % x0, prob % n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = ONE
prob % Delta0 = HALF
prob % calcfc => calcfc_fletcheq2

n = prob % n
call safealloc(prob % lb, n)
prob % lb = -HUGENUM
call safealloc(prob % ub, n)
prob % ub = HUGENUM
call safealloc(prob % Aeq, n, 0_IK)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, n, 0_IK)
call safealloc(prob % bineq, 0_IK)
end subroutine construct_fletcheq2


subroutine calcfc_fletcheq2(x, f, constr)
! Test problem 7 (Equation (14.2.2) in Fletcher's book) in Powell's original COBYLA package.
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: debug_mod, only : assert
implicit none

character(len=*), parameter :: srname = 'CALCFC_FLETCHEQ2'
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: constr(:)
real(RP), intent(out) :: f

call assert(size(x) == 3 .and. size(constr) == 3, 'SIZE(X) == 3, SIZE(CONSTR) == 3', srname)

f = x(3)
constr(1) = 5.0_RP * x(1) - x(2) + x(3)
constr(2) = x(3) - x(1)**2 - x(2)**2 - 4.0 * x(2)
constr(3) = x(3) - 5.0 * x(1) - x(2)
end subroutine calcfc_fletcheq2
