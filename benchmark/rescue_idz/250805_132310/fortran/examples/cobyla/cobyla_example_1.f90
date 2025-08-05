!--------------------------------------------------------------------------------------------------!
! This is an example to illustrate the usage of the solver.
!
! The objective function is trivial. This is intentional, as the focus is how to use the API.
!--------------------------------------------------------------------------------------------------!


!------------------------- THE MODULE THAT IMPLEMENTS CALCFC, CALLBACK_FCN ------------------------!
module calcfc_mod

implicit none
private
public :: IK, RP, calcfc, callback_fcn
integer, parameter :: RP = kind(0.0D0)
integer, parameter :: IK = kind(0)
! N.B.: We assume that PRIMA_REAL_PRECISION = 64 (double precision) and PRIMA_INTEGER_KIND = 0
! (default kind). Revise RP and IK if this is not the case.

contains

! Objective function
! This one is slightly different from the other example_1 files since COBYLA
! requires a different signature for the objective function.
subroutine calcfc(x, f, constr)
implicit none

! Inputs
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)

f = (x(1) - 5.0_RP)**2 + (x(2) - 4.0_RP)**2

! We add a constraint we know will be active in order to demonstrate usage
! The constraint is x(1)**2 - 9 <= 0, meaning |x1| <= 3.
constr(1) = x(1)**2 - 9.0_RP

end subroutine calcfc

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

write (*, '("Best point so far: x = [", F6.4, ", ", F6.4, "], f = ", F6.3, ", cstrv = ", F6.3, &
    & ", nlconstr = [", F6.3, "], nf = ", I0, ", tr = ", I0, "")') x(1), x(2), f, cstrv, nlconstr(1), nf, tr

terminate = .false.

end subroutine callback_fcn

end module calcfc_mod


!---------------------------------------- THE MAIN PROGRAM ----------------------------------------!
program cobyla_exmp

! The following line makes the solver available.
use cobyla_mod, only : cobyla

! The following line specifies which module provides CALCFC and CALLBACK_FCN.
use calcfc_mod, only : RP, IK, calcfc, callback_fcn

implicit none

integer, parameter :: n = 2
integer :: nf, info
real(RP) :: f, x(n), x0(n), cstrv

! Define the starting point.
x0 = 0.0_RP

! The following lines illustrates how to call the solver.
x = x0
call cobyla(calcfc, 1_IK, x, f, cstrv)  ! This call will not print anything.

! In addition to the compulsory arguments, the following illustration specifies also RHOBEG and
! IPRINT, which are optional. All the unspecified optional arguments (RHOEND, MAXFUN, etc.) will
! take their default values coded in the solver.
x = x0
call cobyla(calcfc, 1_IK, x, f, cstrv, rhobeg=1.0_RP, iprint=1_IK, nf=nf, info=info, callback_fcn=callback_fcn)

end program cobyla_exmp
