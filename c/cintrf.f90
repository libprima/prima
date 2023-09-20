module cintrf_mod
!--------------------------------------------------------------------------------------------------!
! cintrf_mod provides interoperability with C objective functions
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: COBJ, COBJCON, evalcobj, evalcobjcon


abstract interface

    subroutine COBJ(x, f) bind(c)
    use, intrinsic :: iso_c_binding, only : C_DOUBLE
    implicit none
    ! We cannot use assumed-shape arrays for C interoperability
    real(C_DOUBLE), intent(in) :: x(*)
    real(C_DOUBLE), intent(out) :: f
    end subroutine COBJ

    subroutine COBJCON(x, f, terminate, constr) bind(c)
    use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_BOOL
    implicit none
    ! We cannot use assumed-shape arrays for C interoperability
    real(C_DOUBLE), intent(in) :: x(*)
    real(C_DOUBLE), intent(out) :: f
    logical(C_BOOL), intent(out) :: terminate
    real(C_DOUBLE), intent(out) :: constr(*)
    end subroutine COBJCON

end interface


contains


subroutine evalcobj(cobj_ptr, x, f)
use, non_intrinsic :: consts_mod, only : RP
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_FUNPTR, C_F_PROCPOINTER
implicit none
type(C_FUNPTR), intent(in) :: cobj_ptr
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f

! Local variables
procedure(COBJ), pointer :: obj_ptr
real(C_DOUBLE) :: x_loc(size(x))
real(C_DOUBLE) :: f_loc

! Read the inputs and convert them to the types specified in COBJ
x_loc = real(x, kind(x_loc))
call C_F_PROCPOINTER(cobj_ptr, obj_ptr)

! Call the C objective function
call obj_ptr(x_loc, f_loc)

! Write the output
f = real(f_loc, kind(f))

end subroutine evalcobj


subroutine evalcobjcon(cobjcon_ptr, x, f, terminate, constr)
use, non_intrinsic :: consts_mod, only : RP
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_FUNPTR, C_F_PROCPOINTER, C_BOOL
implicit none
type(C_FUNPTR), intent(in) :: cobjcon_ptr
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
logical, intent(out) :: terminate
real(RP), intent(out) :: constr(:)

! Local variables
procedure(COBJCON), pointer :: objcon_ptr
real(C_DOUBLE) :: x_loc(size(x))
real(C_DOUBLE) :: f_loc
logical(C_BOOL) :: terminate_loc = .false.
real(C_DOUBLE) :: constr_loc(size(constr))

! Read the inputs and convert them to the types specified in COBJCON
x_loc = real(x, kind(x_loc))
call C_F_PROCPOINTER(cobjcon_ptr, objcon_ptr)

! Call the C objective function
call objcon_ptr(x_loc, f_loc, terminate_loc, constr_loc)

! Write the output
f = real(f_loc, kind(f))
terminate = logical(terminate_loc, kind(terminate))
constr = real(constr_loc, kind(constr))

end subroutine evalcobjcon

end module cintrf_mod
