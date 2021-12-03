module prob_mod
!--------------------------------------------------------------------------------------------------!
! This module implements the following testing problems.
! chebyqad
! chrosen
! hexagon
! trigsabs
! trigssqs
! vardim
!--------------------------------------------------------------------------------------------------!
implicit none
private
public :: PNLEN, probname
public :: getx0, getdelta0, getcsize
public :: calfun, calcfc

integer, parameter :: SEED_DFT = 42
integer, parameter :: PNLEN = 64
character(len=PNLEN) :: probname


contains


function getx0(n) result(x)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, TEN, TENTH, PI
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: noise_mod, only : noisyx
use, non_intrinsic :: param_mod, only : XNOISE_DFT
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
use, non_intrinsic :: string_mod, only : lower, trimstr

implicit none

integer(IK), intent(in) :: n
real(RP) :: x(n)

character(len=*), parameter :: srname = 'GETX0'
integer(IK) :: i
integer, allocatable :: seedsav(:)
real(RP) :: ind(n)
real(RP) :: theta(n)
real(RP) :: xstar(n)
real(RP) :: ystar(n)

ind = real([(i, i=1, n)], RP)

x = ZERO

select case (lower(trimstr(probname)))
case ('chebyqad')
    x = ind / real([(n + 1, i=1, n)], RP)
case ('chrosen')
    x = -ONE
case ('hexagon')
    x = ONE
case ('trigsabs')
    call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
    call setseed(SEED_DFT)  ! Set the random seed by SETSEED(SEED_DFT).
    xstar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{x}^* in the NEWUOA paper.
    ystar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{y}^* in the NEWUOA paper.
    x = xstar + TENTH * ystar
    call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
    deallocate (seedsav)
case ('trigssqs')
    call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
    call setseed(SEED_DFT)  ! Set the random seed by SETSEED(SEED_DFT).
    xstar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{x}^* in the NEWUOA paper.
    ystar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{y}^* in the NEWUOA paper.
    theta = TEN**(-rand(n))
    x = (xstar + TENTH * ystar) / theta
    call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
    deallocate (seedsav)
case ('vardim')
    x = ONE - ind / real(n, RP)
case default
    call errstop(srname, 'Unknown problem: '//trimstr(probname))
end select

x = noisyx(x, noise_level=XNOISE_DFT)

end function getx0


function getdelta0(n) result(delta)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, HALF, TENTH
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: noise_mod, only : noisy
use, non_intrinsic :: param_mod, only : NOISE_DFT
use, non_intrinsic :: string_mod, only : lower, trimstr

implicit none

integer(IK), intent(in) :: n
real(RP) :: delta

character(len=*), parameter :: srname = 'GETDELTA0'

delta = ONE
select case (lower(trimstr(probname)))
case ('chebyqad')
    delta = 0.2_RP / real(n + 1, RP)
case ('chrosen')
    delta = HALF
case ('hexagon')
    delta = HALF
case ('trigsabs')
    delta = TENTH
case ('trigssqs')
    delta = TENTH
case ('vardim')
    delta = ONE / real(2 * n, RP)
case default
    call errstop(srname, 'Unknown problem: '//trimstr(probname))
end select
delta = noisy(delta, noise_level=NOISE_DFT)

end function getdelta0


function getcsize(n) result(m)
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: string_mod, only : lower, trimstr
implicit none

character(len=*), parameter :: srname = 'GETCSIZE'
integer(IK) :: m
integer(IK), intent(in) :: n

select case (lower(trimstr(probname)))
case ('hexagon')
    m = 14
case default
    call errstop(srname, 'Unknown problem: '//trimstr(probname))
end select
end function getcsize


subroutine calfun(x, f)
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: noise_mod, only : noisyfun
use, non_intrinsic :: param_mod, only : FNOISE_DFT
use, non_intrinsic :: string_mod, only : lower, trimstr

implicit none

character(len=*), parameter :: srname = 'CALFUN'
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

select case (lower(trimstr(probname)))
case ('chebyqad')
    f = chebyqad(x)
case ('chrosen')
    f = chrosen(x)
case ('trigsabs')
    f = trigsabs(x)
case ('trigssqs')
    f = trigssqs(x)
case ('vardim')
    f = vardim(x)
case default
    call errstop(srname, 'Unknown problem: '//trimstr(probname))
end select

f = noisyfun(x, f, noise_level=FNOISE_DFT, noise_type='gaussian')

end subroutine calfun


subroutine calcfc(x, f, con)
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: noise_mod, only : noisyfun
use, non_intrinsic :: param_mod, only : FNOISE_DFT
use, non_intrinsic :: string_mod, only : lower, trimstr

implicit none

character(len=*), parameter :: srname = 'CALCFC'
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: con(:)

select case (lower(trimstr(probname)))
case ('hexagon')
    call hexagon(x, f, con)
case default
    call errstop(srname, 'Unknown problem: '//trimstr(probname))
end select

f = noisyfun(x, f, noise_level=FNOISE_DFT, noise_type='gaussian')
con = noisyfun(x, f, noise_level=FNOISE_DFT, noise_type='gaussian')
end subroutine calcfc


pure function chebyqad(x) result(f)
! Chebyquad function (Fletcher, 1965)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO
implicit none

real(RP), intent(in) :: x(:)
real(RP) :: f

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
    if (mod(i, 2_IK) /= 0) then
        tmp = tmp + ONE / real(i * i - 2 * i, RP)
    end if
    f = f + tmp**2
end do

end function chebyqad


pure function chrosen(x) result(f)
! Chained Rosenbrock function.
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
implicit none

real(RP), intent(in) :: x(:)
real(RP) :: f

integer(IK) :: n
real(RP), parameter :: alpha = 1.0E2_RP

n = int(size(x), kind(n))
f = sum((x(1:n - 1) - ONE)**2 + alpha * (x(2:n) - x(1:n - 1)**2)**2); 
end function chrosen


subroutine hexagon(x, f, con)
! Test problem 10 (Hexagon area) in Powell's original COBYLA package.
use, non_intrinsic :: consts_mod, only : RP, ONE, HALF
use, non_intrinsic :: debug_mod, only : assert
implicit none

character(len=*), parameter :: srname = 'HEXAGON'
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: con(:)
real(RP), intent(out) :: f

call assert(size(x) == 9 .and. size(con) == 14, 'SIZE(X) == 9, SIZE(CON) == 14', srname)

f = -HALF * (x(1) * x(4) - x(2) * x(3) + x(3) * x(9) - x(5) * x(9) + x(5) * x(8) - x(6) * x(7))
con(1) = ONE - x(3)**2 - x(4)**2
con(2) = ONE - x(9)**2
con(3) = ONE - x(5)**2 - x(6)**2
con(4) = ONE - x(1)**2 - (x(2) - x(9))**2
con(5) = ONE - (x(1) - x(5))**2 - (x(2) - x(6))**2
con(6) = ONE - (x(1) - x(7))**2 - (x(2) - x(8))**2
con(7) = ONE - (x(3) - x(5))**2 - (x(4) - x(6))**2
con(8) = ONE - (x(3) - x(7))**2 - (x(4) - x(8))**2
con(9) = ONE - x(7)**2 - (x(8) - x(9))**2
con(10) = x(1) * x(4) - x(2) * x(3)
con(11) = x(3) * x(9)
con(12) = -x(5) * x(9)
con(13) = x(5) * x(8) - x(6) * x(7)
con(14) = x(9)
end subroutine hexagon


function trigsabs(x) result(f)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, PI
use, non_intrinsic :: linalg_mod, only : matprod
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
implicit none

real(RP), intent(in) :: x(:)
real(RP) :: f

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
end function trigsabs


function trigssqs(x) result(f)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, TWO, TEN, PI
use, non_intrinsic :: linalg_mod, only : matprod
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
implicit none

real(RP), intent(in) :: x(:)
real(RP) :: f

integer(IK) :: n
integer, allocatable :: seedsav(:)
real(RP) :: C(2 * size(x), size(x))
real(RP) :: S(2 * size(x), size(x))
real(RP) :: theta(size(x))
real(RP) :: xstar(size(x))

n = int(size(x), kind(n))

call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
call setseed(SEED_DFT)  ! Set the random seed by SETSEED(SEED_DFT).
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
end function trigssqs


pure function vardim(x) result(f)
use, non_intrinsic :: consts_mod, only : RP, IK, ONE
implicit none

real(RP), intent(in) :: x(:)
real(RP) :: f

real(RP) :: ind(size(x))
integer(IK) :: i

ind = real([(i, i=1, int(size(x), IK))], RP)
f = sum((x - ONE)**2) + sum(ind * (x - ONE))**2 + sum(ind * (x - ONE))**4
end function vardim


end module prob_mod
