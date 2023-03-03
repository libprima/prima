subroutine construct_vardim(prob, n)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, REALMAX
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
integer(IK) :: i

! Code shared by all unconstrained problems.
prob % probtype = 'u'
prob % m = 0
prob % n = n
call safealloc(prob % lb, n)
prob % lb = -REALMAX
call safealloc(prob % ub, n)
prob % ub = REALMAX
call safealloc(prob % Aeq, n, 0_IK)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, n, 0_IK)
call safealloc(prob % bineq, 0_IK)

! Problem-specific code
prob % probname = 'vardim'
call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = ONE - real([(i, i=1, n)]) / real(n, RP)
prob % Delta0 = ONE / real(2 * n, RP)
prob % calfun => calfun_vardim
prob % calcfc => calcfc_vardim
end subroutine construct_vardim


subroutine calcfc_vardim(x, f, constr)
use, non_intrinsic :: consts_mod, only : RP, ZERO
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
call calfun_vardim(x, f)
constr = ZERO  ! Without this line, compilers may complain that CONSTR is not set.
end subroutine calcfc_vardim


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
