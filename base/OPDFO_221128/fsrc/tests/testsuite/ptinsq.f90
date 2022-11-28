! Test problem in Powell's original BOBYQA package. Powell's description is as follows.
!
! Test problem for BOBYQA, the objective function being the sum of the reciprocals of all pairwise
! distances between the points P_I, I=1,2,...,M in two dimensions, where M=N/2 and where the
! components of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables defines the M points
! P_I. The initial X gives equally spaced points on a circle. Four different choices of the pairs
! (N,NPT) are tried, namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local minimum
! that is not global occurs in both the N=10 cases. The details of the results are highly sensitive
! to computer rounding errors. The choice IPRINT=2 provides the current X and optimal F so far
! whenever RHO is reduced. The bound constraints of the problem require every component of X to be
! in the interval [-1,1].
!
! In Powell's BOBYQA paper, this problem is called "points in square".

subroutine construct_ptinsq(prob, n)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, PI
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Inputs
integer(IK), intent(in) :: n

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
character(len=*), parameter :: srname = 'CONSTRUCT_PTINSQ'
integer(IK) :: j
real(RP) :: angle

call assert(modulo(n, 2_IK) == 0, 'N is even', srname)

prob % probname = 'ptinsq'
prob % probtype = 'b'
prob % n = n

call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
do j = 1_IK, n / 2_IK
    angle = real(j, RP) * TWO * PI / real(n / 2_IK, RP)
    prob % x0(2 * j - 1) = cos(angle)
    prob % x0(2 * j) = sin(angle)
end do

prob % Delta0 = ONE
prob % calfun => calfun_ptinsq

call safealloc(prob % lb, n)
prob % lb = -ONE
call safealloc(prob % ub, n)
prob % ub = ONE

call safealloc(prob % Aeq, n, 0_IK)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, n, 0_IK)
call safealloc(prob % bineq, 0_IK)

end subroutine construct_ptinsq


subroutine calfun_ptinsq(x, f)
use, non_intrinsic :: consts_mod, only : IK, RP
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f

! Local variables
character(len=*), parameter :: srname = 'CALFUN_PTINSQ'
integer(IK) :: i
integer(IK) :: j
integer(IK) :: n
real(RP) :: temp

n = int(size(x), kind(n))

call assert(modulo(n, 2_IK) == 0, 'N is even', srname)

f = 0.0_RP
do i = 4_IK, n, 2_IK
    do j = 2_IK, i - 2_IK, 2_IK
        temp = (x(i - 1) - x(j - 1))**2 + (x(i) - x(j))**2
        temp = max(temp, 1.0E-6_RP)
        f = f + 1.0_RP / sqrt(temp)
    end do
end do

end subroutine calfun_ptinsq
