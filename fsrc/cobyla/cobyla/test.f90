!        This is file : test
! Author= zaikunzhang
! Started at: 29.11.2021
! Last Modified: Thursday, December 02, 2021 AM10:29:53

program test
!use, non_intrinsic :: consts_mod, only : RP, IK
!use, non_intrinsic :: pintrf_mod, only : FUNCON
implicit none

!procedure(FUNCON) :: calcfc
integer, parameter :: IK = kind(1)
integer, parameter :: RP = kind(1.0)
real(RP) :: x(9)
real(RP) :: constr(14)
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP) :: f
real(RP) :: ftarget
real(RP) :: cstrv
real(RP) :: w(3000)
integer(IK) :: iact(51)

real(RP) :: tempa, tempb, tempc, tempd, xopt(size(x))

integer(IK) :: n
integer(IK) :: m
integer(IK) :: iprint
integer(IK) :: maxfun
integer(IK) :: info

iprint = 1
n = 9
m = 14
maxfun = 5000
ftarget = -huge(0.0_RP)
rhobeg = 0.5_RP
rhoend = 0.001_RP
x = 1.0_RP


call COBYLA(N, M, X, RHOBEG, RHOEND, IPRINT, MAXFUN, W, IACT)
write (*, *) info

tempa = x(1) + x(3) + x(5) + x(7)
tempb = x(2) + x(4) + x(6) + x(8)
tempc = 0.5 / sqrt(tempa * tempa + tempb * tempb)
tempd = tempc * sqrt(3.0)
xopt(1) = tempd * tempa + tempc * tempb
xopt(2) = tempd * tempb - tempc * tempa
xopt(3) = tempd * tempa - tempc * tempb
xopt(4) = tempd * tempb + tempc * tempa
xopt(5:8) = xopt(1:4)


write (*, *) cstrv, x
write (*, *) xopt
write (*, *) sum((x - xopt)**2)

!x = 0.0_RP
!rhoend = 0.0001_RP
!!call COBYLA(n, m, x, rhobeg, rhoend, iprint, maxfun, w, iact, f, info, ftarget, cstrv, constr)
!call COBYLA(N, M, X, RHOBEG, RHOEND, IPRINT, MAXFUN, W, IACT)

!tempa = x(1) + x(3) + x(5) + x(7)
!tempb = x(2) + x(4) + x(6) + x(8)
!tempc = 0.5 / sqrt(tempa * tempa + tempb * tempb)
!tempd = tempc * sqrt(3.0)
!xopt(1) = tempd * tempa + tempc * tempb
!xopt(2) = tempd * tempb - tempc * tempa
!xopt(3) = tempd * tempa - tempc * tempb
!xopt(4) = tempd * tempb + tempc * tempa
!xopt(5:8) = xopt(1:4)

!write (*, *) cstrv, x
!write (*, *) xopt
!write (*, *) sum((x - xopt)**2)

end program test
