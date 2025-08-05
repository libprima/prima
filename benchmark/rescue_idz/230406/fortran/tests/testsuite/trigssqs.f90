subroutine construct_trigssqs(prob, n)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, TEN, TENTH, PI, REALMAX
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
integer, allocatable :: seedsav(:)
real(RP) :: theta(n)
real(RP) :: xstar(n)
real(RP) :: ystar(n)

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
prob % probname = 'trigssqs'

call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
call setseed(RANDSEED_DFT)  ! Set the random seed by SETSEED(RANDSEED_DFT).
xstar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{x}^* in the NEWUOA paper.
ystar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{y}^* in the NEWUOA paper.
theta = TEN**(-rand(n))
call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
deallocate (seedsav)
call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = (xstar + TENTH * ystar) / theta

prob % Delta0 = TENTH
prob % calfun => calfun_trigssqs
prob % calcfc => calcfc_trigssqs
end subroutine construct_trigssqs


subroutine calcfc_trigssqs(x, f, constr)
use, non_intrinsic :: consts_mod, only : RP, ZERO
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
call calfun_trigssqs(x, f)
constr = ZERO  ! Without this line, compilers may complain that CONSTR is not set.
end subroutine calcfc_trigssqs


subroutine calfun_trigssqs(x, f)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, TEN, PI
use, non_intrinsic :: linalg_mod, only : matprod
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

integer(IK) :: n
integer, allocatable :: seedsav(:)
real(RP) :: C(2 * size(x), size(x))
real(RP) :: S(2 * size(x), size(x))
real(RP) :: theta(size(x))
real(RP) :: xstar(size(x))

n = int(size(x), kind(n))

call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
call setseed(RANDSEED_DFT)  ! Set the random seed by SETSEED(RANDSEED_DFT).
C = 1.0E2_RP * (TWO * rand(2_IK * n, n) - ONE)
S = 1.0E2_RP * (TWO * rand(2_IK * n, n) - ONE)
xstar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{x}^* in the NEWUOA paper.
! In the NEWUOA paper, THETA is sampled from the logarithmic distribution on [0.1, 1], yet it is not
! clear what Powell meant by "logarithmic distribution". The commonly known logarithmic distribution
! is a discrete distribution on the positive integers.
theta = TEN**(-rand(n))
call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
deallocate (seedsav)

f = sum((matprod(C, cos(xstar) - cos(theta * x)) + matprod(S, sin(xstar) - sin(theta * x)))**2)
end subroutine calfun_trigssqs
