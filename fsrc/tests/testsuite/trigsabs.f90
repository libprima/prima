subroutine construct_trigsabs(prob, n)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, TENTH, PI
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
type(problem_t), intent(out) :: prob

! Local variables
integer, allocatable :: seedsav(:)
real(RP) :: xstar(n)
real(RP) :: ystar(n)

prob % probname = 'trigsabs'
prob % probtype = 'u'
prob % n = n

call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
call setseed(SEED_DFT)  ! Set the random seed by SETSEED(SEED_DFT).
xstar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{x}^* in the NEWUOA paper.
ystar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{y}^* in the NEWUOA paper.
call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
deallocate (seedsav)
call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = xstar + TENTH * ystar

prob % Delta0 = TENTH
prob % calfun => calfun_trigsabs
end subroutine construct_trigsabs


subroutine calfun_trigsabs(x, f)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, PI
use, non_intrinsic :: linalg_mod, only : matprod
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

integer(IK) :: n
integer, allocatable :: seedsav(:)
real(RP) :: C(2 * size(x), size(x))
real(RP) :: S(2 * size(x), size(x))
real(RP) :: xstar(size(x))

n = int(size(x), kind(n))

call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
call setseed(SEED_DFT)  ! Set the random seed by SETSEED(SEED_DFT).
C = 1.0E2_RP * (TWO * rand(2_IK * n, n) - ONE)
S = 1.0E2_RP * (TWO * rand(2_IK * n, n) - ONE)
xstar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{x}^* in the NEWUOA paper.
call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
deallocate (seedsav)

f = sum(abs(matprod(C, cos(xstar) - cos(x)) + matprod(S, sin(xstar) - sin(x))))
end subroutine calfun_trigsabs
