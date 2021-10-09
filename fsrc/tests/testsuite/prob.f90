module prob_mod
implicit none
private
public :: PNLEN, probname, calfun, getx0, getdelta0

integer, parameter :: SEED_DFT = 42
integer, parameter :: PNLEN = 64
character(len=PNLEN) :: probname


contains


function getx0(n) result(x)

use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, TEN, TENTH, PI
use, non_intrinsic :: noise_mod, only : noisyx
use, non_intrinsic :: rand_mod, only : getseed, setseed, rand
use, non_intrinsic :: string_mod, only : lower, trimstr

implicit none

integer(IK), intent(in) :: n
real(RP) :: x(n)

integer(IK) :: i
integer, allocatable :: seedsav(:)
real(RP) :: ind(n)
real(RP) :: theta(n)
real(RP) :: xstar(n)
real(RP) :: ystar(n)

ind = real([(i, i=1, n)], RP)

x = ZERO

if (lower(trimstr(probname)) == 'chebyqad') then
    x = ind / real([(n + 1, i=1, n)], RP)
elseif (lower(trimstr(probname)) == 'chrosen') then
    x = -ONE
elseif (lower(trimstr(probname)) == 'trigsabs') then
    call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
    call setseed(SEED_DFT)  ! Set the random seed by SETSEED(SEED_DFT).
    xstar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{x}^* in the NEWUOA paper.
    ystar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{y}^* in the NEWUOA paper.
    x = xstar + TENTH * ystar
    call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
    deallocate (seedsav)
elseif (lower(trimstr(probname)) == 'trigssqs') then
    call getseed(seedsav)  ! Backup the current random seed in SEEDSAV.
    call setseed(SEED_DFT)  ! Set the random seed by SETSEED(SEED_DFT).
    xstar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{x}^* in the NEWUOA paper.
    ystar = PI * (TWO * rand(n) - ONE)  ! This is the \hat{y}^* in the NEWUOA paper.
    theta = TEN**(-rand(n))
    x = (xstar + TENTH * ystar) / theta
    call setseed(seedsav)  ! Recover the random seed by SEEDSAV.
    deallocate (seedsav)
elseif (lower(trimstr(probname)) == 'vardim') then
    x = ONE - ind / real(n, RP)
end if

x = noisyx(x)

end function getx0


function getdelta0(n) result(delta)

use, non_intrinsic :: consts_mod, only : RP, IK, ONE, HALF, TENTH
use, non_intrinsic :: noise_mod, only : noisy
use, non_intrinsic :: string_mod, only : lower, trimstr

implicit none

integer(IK), intent(in) :: n
real(RP) :: delta

delta = ONE
if (lower(trimstr(probname)) == 'chebyqad') then
    delta = 0.2_RP / real(n + 1, RP)
elseif (lower(trimstr(probname)) == 'chrosen') then
    delta = HALF
elseif (lower(trimstr(probname)) == 'trigsabs') then
    delta = TENTH
elseif (lower(trimstr(probname)) == 'trigssqs') then
    delta = TENTH
elseif (lower(trimstr(probname)) == 'vardim') then
    delta = ONE / real(2 * n, RP)
end if
delta = noisy(delta)

end function getdelta0


subroutine calfun(x, f)

use, non_intrinsic :: consts_mod, only : RP, ZERO
use, non_intrinsic :: noise_mod, only : noisyfun
use, non_intrinsic :: string_mod, only : lower, trimstr

implicit none

real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

f = ZERO

if (lower(trimstr(probname)) == 'chebyqad') then
    f = chebyqad(x)
elseif (lower(trimstr(probname)) == 'chrosen') then
    f = chrosen(x)
elseif (lower(trimstr(probname)) == 'trigsabs') then
    f = trigsabs(x)
elseif (lower(trimstr(probname)) == 'trigssqs') then
    f = trigssqs(x)
elseif (lower(trimstr(probname)) == 'vardim') then
    f = vardim(x)
end if

f = noisyfun(x, f)

end subroutine calfun


pure function chebyqad(x) result(f)
! Chebyquad function (Fletcher, 1965)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO
implicit none

real(RP), intent(in) :: x(:)
real(RP) :: f

integer(IK) :: i, n
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
