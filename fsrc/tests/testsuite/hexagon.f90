subroutine construct_hexagon(prob)
use, non_intrinsic :: consts_mod, only : ONE, HALF, HUGENUM
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Outputs
type(problem_t), intent(out) :: prob

! Local variables
integer(IK) :: n

prob % probname = 'hexagon'
prob % probtype = 'n'
prob % m = 14
prob % n = 9
call safealloc(prob % x0, prob % n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = ONE
prob % Delta0 = HALF
prob % calcfc => calcfc_hexagon

n = prob % n
call safealloc(prob % lb, n)
prob % lb = -HUGENUM
call safealloc(prob % ub, n)
prob % ub = HUGENUM
call safealloc(prob % Aeq, 0_IK, n)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, 0_IK, n)
call safealloc(prob % bineq, 0_IK)
end subroutine construct_hexagon


subroutine calcfc_hexagon(x, f, con)
! Test problem 10 (Hexagon area) in Powell's original COBYLA package.
use, non_intrinsic :: consts_mod, only : RP, ONE, HALF
use, non_intrinsic :: debug_mod, only : assert
implicit none

character(len=*), parameter :: srname = 'CALCFC_HEXAGON'
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: con(:)
real(RP), intent(out) :: f

call assert(size(x) == 9 .and. size(con) == 14, 'SIZE(X) == 9, SIZE(CON) == 14', srname)

f = -HALF * (x(1) * x(4) - x(2) * x(3) + x(3) * x(9) - x(5) * x(9) + x(5) * x(8) - x(6) * x(7))
con(1) = ONE - x(3)**2 - x(4)**2
con(2) = ONE - x(9)**2
con(3) = ONE - x(5)**2 - x(6)**2
con(4) = ONE - x(1)**2 - (x(2) - x(9))**2
con(5) = ONE - (x(1) - x(5))**2 - (x(2) - x(6))**2
con(6) = ONE - (x(1) - x(7))**2 - (x(2) - x(8))**2
con(7) = ONE - (x(3) - x(5))**2 - (x(4) - x(6))**2
con(8) = ONE - (x(3) - x(7))**2 - (x(4) - x(8))**2
con(9) = ONE - x(7)**2 - (x(8) - x(9))**2
con(10) = x(1) * x(4) - x(2) * x(3)
con(11) = x(3) * x(9)
con(12) = -x(5) * x(9)
con(13) = x(5) * x(8) - x(6) * x(7)
con(14) = x(9)
end subroutine calcfc_hexagon
