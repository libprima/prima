subroutine construct_vardim(prob, n)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
type(problem_t), intent(out) :: prob

! Local variables
integer(IK) :: i

prob % probname = 'vardim'
prob % probtype = 'u'
prob % n = n
call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = ONE - real([(i, i=1, n)]) / real(n, RP)
prob % Delta0 = ONE / real(2 * n, RP)
prob % calfun => calfun_vardim
end subroutine construct_vardim


subroutine calfun_vardim(x, f)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

real(RP) :: ind(size(x))
integer(IK) :: i

ind = real([(i, i=1, int(size(x), IK))], RP)
f = sum((x - ONE)**2) + sum(ind * (x - ONE))**2 + sum(ind * (x - ONE))**4
end subroutine calfun_vardim
