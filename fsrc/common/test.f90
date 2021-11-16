!        This is file : test
! Author= zaikunzhang
! Started at: 10.11.2021
! Last Modified: Tuesday, November 16, 2021 PM12:07:01

program test
use consts_mod
use linalg_mod
!use rand_mod
implicit none

integer, parameter :: n = 5000
real(RP) :: G(2, 2)
real(RP), parameter :: a = sqrt(huge(0.0_RP))
!real(RP), parameter :: a = sqrt(tiny(0.0_RP))
!real(RP), parameter :: a = sqrt(epsilon(0.0_RP))
real(RP), parameter :: x(2) = [a / PI, a / (PI * 1.0E2_RP)]
real(RP), parameter :: y(2) = [x(2), x(1)]
real(RP) :: r

r = sqrt(sum(x**2))
G = planerot(x)
write (*, *) '-->', maxval(abs(matprod(transpose(G), G) - eye(2))), norm(matprod(G, x) - [r, ZERO]) / max(r, EPS)
G = planerot(y)
write (*, *) '-->', maxval(abs(matprod(transpose(G), G) - eye(2))), norm(matprod(G, y) - [r, ZERO]) / max(r, EPS)

end program test
