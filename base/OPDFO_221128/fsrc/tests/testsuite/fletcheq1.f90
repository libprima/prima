subroutine construct_fletcheq1(prob)
use, non_intrinsic :: consts_mod, only : ONE, HALF, HUGENUM
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
integer(IK) :: n

prob % probname = 'fletcheq1'
prob % probtype = 'n'
prob % m = 2
prob % n = 2
call safealloc(prob % x0, prob % n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = ONE
prob % Delta0 = HALF
prob % calcfc => calcfc_fletcheq1

n = prob % n
call safealloc(prob % lb, n)
prob % lb = -HUGENUM
call safealloc(prob % ub, n)
prob % ub = HUGENUM
call safealloc(prob % Aeq, n, 0_IK)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, n, 0_IK)
call safealloc(prob % bineq, 0_IK)
end subroutine construct_fletcheq1


subroutine calcfc_fletcheq1(x, f, constr)
! Test problem 6 (Equation (9.1.15) in Fletcher's book) in Powell's original COBYLA package.
use, non_intrinsic :: consts_mod, only : RP, ONE
use, non_intrinsic :: debug_mod, only : assert
implicit none

character(len=*), parameter :: srname = 'CALCFC_FLETCHEQ1'
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: constr(:)
real(RP), intent(out) :: f

call assert(size(x) == 2 .and. size(constr) == 2, 'SIZE(X) == 2, SIZE(CONSTR) == 2', srname)

f = -x(1) - x(2)
constr(1) = x(2) - x(1)**2
constr(2) = ONE - x(1)**2 - x(2)**2
end subroutine calcfc_fletcheq1
