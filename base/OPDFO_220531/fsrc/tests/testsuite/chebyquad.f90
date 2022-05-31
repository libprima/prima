subroutine construct_chebyquad(prob, n)
use, non_intrinsic :: consts_mod, only : RP, IK, HUGENUM
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
prob % lb = -HUGENUM
call safealloc(prob % ub, n)
prob % ub = HUGENUM
call safealloc(prob % Aeq, n, 0_IK)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, n, 0_IK)
call safealloc(prob % bineq, 0_IK)

! Problem-specific code
prob % probname = 'chebyquad'
call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0 = real([(i, i=1, n)], RP) / real(n + 1, RP)
prob % Delta0 = 0.2_RP / real(n + 1, RP)
prob % calfun => calfun_chebyquad
prob % calcfc => calcfc_chebyquad
end subroutine construct_chebyquad


subroutine calcfc_chebyquad(x, f, constr)
use, non_intrinsic :: consts_mod, only : RP, ZERO
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
call calfun_chebyquad(x, f)
constr = ZERO  ! Without this line, compilers may complain that CONSTR is not set.
end subroutine calcfc_chebyquad


subroutine calfun_chebyquad(x, f)
! Chebyquad function (Fletcher, 1965)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f

! Local variables
integer(IK) :: i
integer(IK) :: n
real(RP) :: y(size(x) + 1, size(x) + 1), tmp

n = int(size(x), kind(n))

y(1:n, 1) = ONE
y(1:n, 2) = TWO * x - ONE
do i = 2, n
    y(1:n, i + 1) = TWO * y(1:n, 2) * y(1:n, i) - y(1:n, i - 1)
end do

f = ZERO
do i = 1, n + 1_IK
    tmp = sum(y(1:n, i)) / real(n, RP)
    if (modulo(i, 2_IK) /= 0) then
        tmp = tmp + ONE / real(i * i - 2 * i, RP)
    end if
    f = f + tmp**2
end do

end subroutine calfun_chebyquad
