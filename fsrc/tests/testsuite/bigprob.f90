!--------------------------------------------------------------------------------------------------!
! This is a "big" problem that is intended to test whether the solvers can encounter stack overflow.
!--------------------------------------------------------------------------------------------------!

subroutine construct_bigprob(prob, n, m)
use, non_intrinsic :: consts_mod, only : IK, ONE, HALF, TEN
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
integer(IK), intent(in), optional :: m
integer(IK), intent(in) :: n

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
integer(IK) :: i, j, m_loc

if (present(m)) then
    m_loc = m
else
    m_loc = int(min(50 * int(n), 10**(range(0_IK) / 2)), IK)
end if

prob % probtype = 'n'
prob % n = n
prob % m = m_loc

call safealloc(prob % lb, n)
prob % lb = -TEN
call safealloc(prob % ub, n)
prob % ub = TEN

call safealloc(prob % Aeq, n, m_loc)
do j = 1, m_loc
    do i = 1, n
        prob % Aeq(i, j) = ONE / real(i + j, RP)
    end do
end do

call safealloc(prob % beq, m_loc)
prob % beq = sin(real([(j, j=1, m_loc)], RP))

call safealloc(prob % Aineq, n, m_loc)
prob % Aineq = prob % Aeq

call safealloc(prob % bineq, m_loc)
prob % bineq = prob % beq + ONE

! Problem-specific code
prob % probname = 'bigprob'
call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = -ONE
prob % Delta0 = HALF
prob % calfun => calfun_bigprob
prob % calcfc => calcfc_bigprob
end subroutine construct_bigprob


subroutine calcfc_bigprob(x, f, constr)
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
integer(IK) :: i, m
call calfun_bigprob(x, f)
m = int(size(constr), kind(m))
constr = sin(f + sum(abs(x)) + real([(i, i=1, m)], RP))
end subroutine calcfc_bigprob


subroutine calfun_bigprob(x, f)
! Chained Rosenbrock function.
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

integer(IK) :: n
real(RP), parameter :: alpha = 1.0_RP

n = int(size(x), kind(n))
f = sum((x(1:n - 1) - ONE)**2 + alpha * (x(2:n) - x(1:n - 1)**2)**2); 
end subroutine calfun_bigprob
