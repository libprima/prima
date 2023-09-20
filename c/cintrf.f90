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

    subroutine COBJ(x, f, data_ptr) bind(c)
    use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_PTR
    implicit none
    ! We cannot use assumed-shape arrays for C interoperability
    real(C_DOUBLE), intent(in) :: x(*)
    real(C_DOUBLE), intent(out) :: f
    type(C_PTR), intent(in), value :: data_ptr
    end subroutine COBJ

    subroutine COBJCON(x, f, constr, data_ptr) bind(c)
    use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_PTR
    implicit none
    ! We cannot use assumed-shape arrays for C interoperability
    real(C_DOUBLE), intent(in) :: x(*)
    real(C_DOUBLE), intent(out) :: f
    real(C_DOUBLE), intent(out) :: constr(*)
    type(C_PTR), intent(in), value :: data_ptr
    end subroutine COBJCON

end interface


contains


subroutine evalcobj(cobj_ptr, data_ptr, x, f)
use, non_intrinsic :: consts_mod, only : RP
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_FUNPTR, C_F_PROCPOINTER, C_PTR
implicit none
type(C_FUNPTR), intent(in) :: cobj_ptr
type(C_PTR), intent(in), value :: data_ptr
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
call obj_ptr(x_loc, f_loc, data_ptr)

! Write the output
f = real(f_loc, kind(f))

end subroutine evalcobj


subroutine evalcobjcon(cobjcon_ptr, data_ptr, x, f, constr)
use, non_intrinsic :: consts_mod, only : RP
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_FUNPTR, C_F_PROCPOINTER, C_PTR
implicit none
type(C_FUNPTR), intent(in) :: cobjcon_ptr
type(C_PTR), intent(in), value :: data_ptr
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)

! Local variables
procedure(COBJCON), pointer :: objcon_ptr
real(C_DOUBLE) :: x_loc(size(x))
real(C_DOUBLE) :: f_loc
real(C_DOUBLE) :: constr_loc(size(constr))

! Read the inputs and convert them to the types specified in COBJCON
x_loc = real(x, kind(x_loc))
call C_F_PROCPOINTER(cobjcon_ptr, objcon_ptr)

! Call the C objective function
call objcon_ptr(x_loc, f_loc, constr_loc, data_ptr)

! Write the output
f = real(f_loc, kind(f))
constr = real(constr_loc, kind(constr))

end subroutine evalcobjcon

end module cintrf_mod
