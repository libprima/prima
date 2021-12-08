subroutine construct_chrosen(prob, n)
use, non_intrinsic :: consts_mod, only : IK, ONE, HALF
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
type(problem_t), intent(out) :: prob

prob % probname = 'chrosen'
prob % probtype = 'u'
prob % n = n
call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = -ONE
prob % Delta0 = HALF
prob % calfun => calfun_chrosen
end subroutine construct_chrosen


subroutine calfun_chrosen(x, f)
! Chained Rosenbrock function.
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

integer(IK) :: n
real(RP), parameter :: alpha = 1.0E2_RP

n = int(size(x), kind(n))
f = sum((x(1:n - 1) - ONE)**2 + alpha * (x(2:n) - x(1:n - 1)**2)**2); 
end subroutine calfun_chrosen
