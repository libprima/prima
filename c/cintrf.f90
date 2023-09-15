module cintrf_mod
!--------------------------------------------------------------------------------------------------!
! cintrf_mod provides interoperability with C objective functions
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: COBJ, COBJCON, evalcobj, evalcobjcon, evalcallback


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

    subroutine CCALLBACK(n, x, f, nf, tr, cstrv, m_nlcon, nlconstr, terminate) bind(c)
    use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_BOOL, C_INT
    implicit none
    integer(C_INT), intent(in), value :: n
    ! We cannot use assumed-shape arrays for C interoperability
    real(C_DOUBLE), intent(in) :: x(*)
    real(C_DOUBLE), intent(in), value :: f
    integer(C_INT), intent(in), value :: nf
    integer(C_INT), intent(in), value :: tr
    real(C_DOUBLE), intent(in), value :: cstrv
    integer(C_INT), intent(in), value :: m_nlcon
    real(C_DOUBLE), intent(in) :: nlconstr(*)
    logical(C_BOOL), intent(out) :: terminate
    end subroutine CCALLBACK


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


subroutine evalcallback(fcb_ptr, x, f, nf, tr, cstrv, nlconstr, terminate)
use, non_intrinsic :: consts_mod, only : RP, IK
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_FUNPTR, C_F_PROCPOINTER, C_BOOL, C_INT
implicit none
type(C_FUNPTR), intent(in) :: fcb_ptr
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: f
integer(IK), intent(in) :: nf
integer(IK), intent(in) :: tr
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: nlconstr(:)
logical, intent(out) :: terminate

! Local variables
procedure(CCALLBACK), pointer :: cb_ptr
integer(C_INT) :: n_loc
real(C_DOUBLE) :: x_loc(size(x))
real(C_DOUBLE) :: f_loc
integer(C_INT) :: nf_loc
integer(C_INT) :: tr_loc
real(C_DOUBLE) :: cstrv_loc
real(C_DOUBLE) :: nlconstr_loc(size(nlconstr))
logical(C_BOOL) :: terminate_loc

! Read the inputs and convert them to the types specified in CCALLBACK
n_loc = size(x)
x_loc = real(x, kind(x_loc))
f_loc = real(f, kind(f_loc))
nf_loc = int(nf, kind(nf_loc))
tr_loc = int(tr, kind(tr_loc))
cstrv_loc = real(cstrv, kind(cstrv_loc))
nlconstr_loc = real(nlconstr, kind(nlconstr_loc))
call C_F_PROCPOINTER(fcb_ptr, cb_ptr)

! Call the C objective function
call cb_ptr(n_loc, x_loc, f_loc, nf_loc, tr_loc, cstrv_loc, size(nlconstr_loc), nlconstr_loc, terminate_loc)

! Write the output
terminate = logical(terminate_loc, kind(terminate))

end subroutine evalcallback


end module cintrf_mod
