subroutine construct_rsnszk(prob)
use, non_intrinsic :: consts_mod, only : ONE, HALF, HUGENUM
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
integer(IK) :: n

prob % probname = 'rsnszk'
prob % probtype = 'n'
prob % m = 3
prob % n = 4
call safealloc(prob % x0, prob % n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = ONE
prob % Delta0 = HALF
prob % calcfc => calcfc_rsnszk

n = prob % n
call safealloc(prob % lb, n)
prob % lb = -HUGENUM
call safealloc(prob % ub, n)
prob % ub = HUGENUM
call safealloc(prob % Aeq, n, 0_IK)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, n, 0_IK)
call safealloc(prob % bineq, 0_IK)
end subroutine construct_rsnszk


subroutine calcfc_rsnszk(x, f, constr)
! Test problem 8 (Rosen-Suzuki) in Powell's original COBYLA package.
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: debug_mod, only : assert
implicit none

character(len=*), parameter :: srname = 'CALCFC_RSNSZK'
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: constr(:)
real(RP), intent(out) :: f

call assert(size(x) == 4 .and. size(constr) == 3, 'SIZE(X) == 4, SIZE(CONSTR) == 3', srname)

f = x(1)**2 + x(2)**2 + 2.0_RP * x(3)**2 + x(4)**2 - 5.0_RP * x(1) - 5.0_RP * x(2) - 21.0_RP * x(3) + 7.0_RP * x(4)
constr(1) = 8.0_RP - x(1)**2 - x(2)**2 - x(3)**2 - x(4)**2 - x(1) + x(2) - x(3) + x(4)
constr(2) = 10.0_RP - x(1)**2 - 2.0_RP * x(2)**2 - x(3)**2 - 2.0_RP * x(4)**2 + x(1) + x(4)
constr(3) = 5.0_RP - 2.0_RP * x(1)**2 - x(2)**2 - x(3)**2 - 2.0_RP * x(1) + x(2) + x(4)
end subroutine calcfc_rsnszk
