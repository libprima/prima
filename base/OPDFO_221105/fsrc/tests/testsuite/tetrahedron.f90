! Test problem in Powell's original LINCOA package. Powell's description is as follows.
!
! Calculate the tetrahedron of least volume that encloses the points (XP(J),YP(J),ZP(J)),
! J=1,2,...,NP. Our method requires the origin to be strictly inside the convex hull of these
! points. There are twelve variables that define the four faces of each tetrahedron that is
! considered. Each face has the form ALPHA*X+BETA*Y+GAMMA*Z=1, the variables X(3K-2), X(3K-1) and
! X(3K) being the values of ALPHA, BETA and GAMMA for the K-th face, K=1,2,3,4. Let the set
! T contain all points in three dimensions that can be reached from the origin without crossing
! a face. Because the volume of T may be infinite, the objective function is the smaller of FMAX and
! the volume of T, where FMAX is set to an upper bound on the final volume initially.  There are
! 4*NP linear constraints on the variables, namely that each of the given points (XP(J),YP(J),ZP(J))
! shall be in T. Let XS = min XP(J), YS = min YP(J), ZS = min ZP(J) and SS = max XP(J)+YP(J)+ZP(J),
! where J runs from 1 to NP. The initial values of the variables are X(1)=1/XS, X(5)=1/YS,
! X(9)=1/ZS, X(2)=X(3)=X(4)=X(6)=X(7) =X(8)=0 and X(10)=X(11)=X(12)=1/SS, which satisfy the linear
! constraints, and which provide the bound FMAX=(SS-XS-YS-ZS)**3/6. Other details of the test
! calculation are given below, including the choice of the data points (XP(J),YP(J),ZP(J)),
! J=1,2,...,NP.

subroutine construct_tetrahedron(prob)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, PI, HUGENUM
use, non_intrinsic :: memory_mod, only : safealloc
implicit none

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
integer(IK) :: i
integer(IK) :: j
integer(IK), parameter :: n = 12_IK
integer(IK), parameter :: np = 50_IK
real(RP) :: Aineq(n, 4 * np)
real(RP) :: ss
real(RP) :: theta(np)
real(RP) :: xp(np)
real(RP) :: xs
real(RP) :: yp(np)
real(RP) :: ys
real(RP) :: zp(np)
real(RP) :: zs

prob % probname = 'tetrahedron'
prob % probtype = 'l'
prob % n = n

! Set the parameters needed for defining X0, Aineq, and bineq.
theta = PI*[(real(j - 1_IK, RP) / real(np - 1_IK, RP), j=1, np)]
xp = cos(theta) * cos(TWO * theta)
yp = sin(theta) * cos(TWO * theta)
zp = sin(TWO * theta)
xp = xp - sum(xp) / real(np, RP)
yp = yp - sum(yp) / real(np, RP)
zp = zp - sum(zp) / real(np, RP)
xs = minval([ZERO, xp])
ys = minval([ZERO, yp])
zs = minval([ZERO, zp])
ss = maxval([ZERO, xp + yp + zp])

call safealloc(prob % x0, n)  ! Not needed if F2003 is fully supported. Needed by Absoft 22.0.
prob % x0(2:8) = ZERO
prob % x0(1) = ONE / xs
prob % x0(5) = ONE / ys
prob % x0(9) = ONE / zs
prob % x0(10:12) = ONE / ss

prob % Delta0 = ONE
prob % calfun => calfun_tetrahedron

call safealloc(prob % lb, n)
prob % lb = -HUGENUM
call safealloc(prob % ub, n)
prob % ub = HUGENUM

call safealloc(prob % Aeq, n, 0_IK)
call safealloc(prob % beq, 0_IK)
call safealloc(prob % Aineq, n, 4_IK * np)
call safealloc(prob % bineq, 4_IK * np)

Aineq = ZERO
do j = 1_IK, np
    do i = 1_IK, 4_IK
        Aineq(3_IK * i - 2_IK, 4_IK * j + i - 4_IK) = xp(j)
        Aineq(3_IK * i - 1_IK, 4_IK * j + i - 4_IK) = yp(j)
        Aineq(3_IK * i, 4_IK * j + i - 4_IK) = zp(j)
    end do
end do
prob % Aineq = Aineq
prob % bineq = ONE

end subroutine construct_tetrahedron


subroutine calfun_tetrahedron(x, f)
use, non_intrinsic :: consts_mod, only : IK, RP, ZERO, TWO, PI
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f

! Local variables
character(len=*), parameter :: srname = 'CALFUN_TETRAHEDRON'
integer(IK) :: j
integer(IK), parameter :: np = 50_IK
real(RP) :: ss
real(RP) :: theta(np)
real(RP) :: v12, v13, v14, v23, v24, v34, del1, del2, del3, del4, temp
real(RP) :: xp(np)
real(RP) :: xs
real(RP) :: yp(np)
real(RP) :: ys
real(RP) :: zp(np)
real(RP) :: zs

call assert(size(x) == 12, 'SIZE(X) == 12', srname)

theta = PI*[(real(j - 1_IK, RP) / real(np - 1_IK, RP), j=1_IK, np)]
xp = cos(theta) * cos(TWO * theta)
yp = sin(theta) * cos(TWO * theta)
zp = sin(TWO * theta)
xp = xp - sum(xp) / real(np, RP)
yp = yp - sum(yp) / real(np, RP)
zp = zp - sum(zp) / real(np, RP)
xs = minval([ZERO, xp])
ys = minval([ZERO, yp])
zs = minval([ZERO, zp])
ss = maxval([ZERO, xp + yp + zp])
f = (ss - xs - ys - zs)**3 / 6.0_RP

v12 = x(1) * x(5) - x(4) * x(2)
v13 = x(1) * x(8) - x(7) * x(2)
v14 = x(1) * x(11) - x(10) * x(2)
v23 = x(4) * x(8) - x(7) * x(5)
v24 = x(4) * x(11) - x(10) * x(5)
v34 = x(7) * x(11) - x(10) * x(8)
del1 = v23 * x(12) - v24 * x(9) + v34 * x(6)
if (del1 <= 0) return
del2 = -v34 * x(3) - v13 * x(12) + v14 * x(9)
if (del2 <= 0) return
del3 = -v14 * x(6) + v24 * x(3) + v12 * x(12)
if (del3 <= 0) return
del4 = -v12 * x(9) + v13 * x(6) - v23 * x(3)
if (del4 <= 0) return
temp = (del1 + del2 + del3 + del4)**3 / (del1 * del2 * del3 * del4)
f = min(temp / 6.0_RP, f)

end subroutine calfun_tetrahedron
