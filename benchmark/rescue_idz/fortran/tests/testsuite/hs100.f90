subroutine construct_hs100(prob)
use, non_intrinsic :: consts_mod, only : ONE, HALF, REALMAX
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
integer(IK) :: n

prob % probname = 'hs100'
prob % probtype = 'n'
prob % m = 4
prob % n = 7
call safealloc(prob % x0, prob % n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = ONE
prob % Delta0 = HALF
prob % calcfc => calcfc_hs100

n = prob % n
call safealloc(prob % lb, n)
prob % lb = -REALMAX
call safealloc(prob % ub, n)
prob % ub = REALMAX
call safealloc(prob % Aeq, n, 0_IK)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, n, 0_IK)
call safealloc(prob % bineq, 0_IK)
end subroutine construct_hs100


subroutine calcfc_hs100(x, f, constr)
! Test problem 9 (HS100) in Powell's original COBYLA package.
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: debug_mod, only : assert
implicit none

character(len=*), parameter :: srname = 'CALCFC_HS100'
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: constr(:)
real(RP), intent(out) :: f

call assert(size(x) == 7 .and. size(constr) == 4, 'SIZE(X) == 7, SIZE(CONSTR) == 4', srname)

f = (x(1) - 10.0_RP)**2 + 5.0_RP * (x(2) - 12.0_RP)**2 + x(3)**4 + 3.0_RP * (x(4) - 11.0_RP)**2 + &
    & 10.0_RP * x(5)**6 + 7.0_RP * x(6)**2 + x(7)**4 - 4.0_RP * x(6) * x(7) - 10.0_RP * x(6) - 8.0_RP * x(7)
constr(1) = 127.0_RP - 2.0_RP * x(1)**2 - 3.0_RP * x(2)**4 - x(3) - 4.0_RP * x(4)**2 - 5.0_RP * x(5)
constr(2) = 282.0_RP - 7.0_RP * x(1) - 3.0_RP * x(2) - 10.0_RP * x(3)**2 - x(4) + x(5)
constr(3) = 196.0_RP - 23.0_RP * x(1) - x(2)**2 - 6.0_RP * x(6)**2 + 8.0_RP * x(7)
constr(4) = -4.0_RP * x(1)**2 - x(2)**2 + 3.0_RP * x(1) * x(2) - 2.0_RP * x(3)**2 - 5.0_RP * x(6) + 11.0_RP * x(7)
end subroutine calcfc_hs100
