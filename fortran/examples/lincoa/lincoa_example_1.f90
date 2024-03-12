!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! The objective function is trivial. This is intentional, as the focus is how to use the API.
!--------------------------------------------------------------------------------------------------!


!------------------------- THE MODULE THAT IMPLEMENTS CALFUN, CALLBACK_FCN ------------------------!
module calfun_mod

implicit none
private
public :: IK, RP, calfun, callback_fcn
integer, parameter :: RP = kind(0.0D0)
integer, parameter :: IK = kind(0)
! N.B.: We assume that PRIMA_REAL_PRECISION = 64 (double precision) and PRIMA_INTEGER_KIND = 0
! (default kind). Revise RP and IK if this is not the case.

contains

! Objective function
subroutine calfun(x, f)
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f

f = (x(1) - 5.0_RP)**2 + (x(2) - 4.0_RP)**2

end subroutine calfun

! Callback function
subroutine callback_fcn(x, f, nf, tr, cstrv, nlconstr, terminate)
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: f
integer(IK), intent(in) :: nf
integer(IK), intent(in) :: tr
real(RP), intent(in), optional :: cstrv
real(RP), intent(in), optional :: nlconstr(:)
logical, intent(out), optional :: terminate

if (.false.) print *, nlconstr   ! Suppress compiler warning about unused variable

write (*, '("Best point so far: x = [", F6.4, ", ", F6.4, "], f = ", F6.3, ", cstrv = ", F6.3, &
&    ", nf = ", I0, ", tr = ", I0, "")') x(1), x(2), f, cstrv, nf, tr

terminate = .false.

end subroutine callback_fcn

end module calfun_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program lincoa_exmp

! The following line makes the solver available.
use lincoa_mod, only : lincoa

! The following line specifies which module provides CALFUN and CALLBACK_FCN.
use calfun_mod, only : RP, IK, calfun, callback_fcn

implicit none

integer, parameter :: n = 2
integer :: nf, info
real(RP) :: f, cstrv, x(n), x0(n), Aineq(1, n), bineq(1)

! Define the starting point.
x0 = 0.0_RP

! Define the constraints. We define constraints that will be active
! in order to demonstrate their usage. The constraint is x1 + x2 < 8
Aineq = 1.0_RP
bineq = 8.0_RP

! The following lines illustrates how to call the solver.
x = x0
call lincoa(calfun, x, f, cstrv, Aineq, bineq)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
x = x0
call lincoa(calfun, x, f, cstrv, Aineq, bineq, rhobeg=1.0_RP, iprint=1_IK, nf=nf, info=info, callback_fcn=callback_fcn)

end program lincoa_exmp
