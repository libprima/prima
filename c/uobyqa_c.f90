module uobyqa_c_mod
!--------------------------------------------------------------------------------------------------!
! uobyqa_c_mod provides uobyqa_c, a simplified interface to cobyla for interoperability with C
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!--------------------------------------------------------------------------------------------------!
implicit none
private
public :: uobyqa_c


contains


subroutine uobyqa_c(cobj_ptr, data_ptr, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, iprint, callback_ptr, info) bind(C)
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_INT, C_FUNPTR, C_PTR, C_ASSOCIATED, C_F_PROCPOINTER
use, non_intrinsic :: cintrf_mod, only : COBJ, CCALLBACK
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: uobyqa_mod, only : uobyqa
implicit none

! Compulsory arguments
type(C_FUNPTR), intent(IN), value :: cobj_ptr
type(C_PTR), intent(in), value :: data_ptr

integer(C_INT), intent(in), value :: n
! We cannot use assumed-shape arrays for C interoperability
real(C_DOUBLE), intent(inout) :: x(n)
real(C_DOUBLE), intent(out) :: f
integer(C_INT), intent(out) :: nf
real(C_DOUBLE), intent(in), value :: rhobeg
real(C_DOUBLE), intent(in), value :: rhoend
real(C_DOUBLE), intent(in), value :: ftarget
integer(C_INT), intent(in), value :: maxfun
integer(C_INT), intent(in), value :: iprint
type(C_FUNPTR), intent(in), value :: callback_ptr
integer(C_INT), intent(out) :: info

! Local variables
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfun_loc
integer(IK) :: nf_loc
real(RP) :: f_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP) :: ftarget_loc
real(RP) :: x_loc(n)
procedure(CCALLBACK), pointer :: cb_ptr
procedure(COBJ), pointer :: obj_ptr

! Read the inputs and convert them to the Fortran side types
x_loc = real(x, kind(x_loc))
rhobeg_loc = real(rhobeg, kind(rhobeg_loc))
rhoend_loc = real(rhoend, kind(rhoend_loc))
ftarget_loc = real(ftarget, kind(ftarget_loc))
maxfun_loc = int(maxfun, kind(maxfun_loc))
iprint_loc = int(iprint, kind(iprint_loc))
call C_F_PROCPOINTER(cobj_ptr, obj_ptr)

! Call the Fortran code
if (C_ASSOCIATED(callback_ptr)) then
    ! If a C callback function is provided, we capture it for use in the closure below
    call C_F_PROCPOINTER(callback_ptr, cb_ptr)
    ! And then we pass the closure to the Fortran code
    call uobyqa(calfun, x_loc, f_loc, nf=nf_loc, rhobeg=rhobeg_loc, rhoend=rhoend_loc, ftarget=ftarget_loc, &
        & maxfun=maxfun_loc, iprint=iprint_loc, callback_fcn=callback_fcn, info=info_loc)
else
    call uobyqa(calfun, x_loc, f_loc, nf=nf_loc, rhobeg=rhobeg_loc, rhoend=rhoend_loc, ftarget=ftarget_loc, &
        & maxfun=maxfun_loc, iprint=iprint_loc, info=info_loc)
end if

! Write the outputs
x = real(x_loc, kind(x))
f = real(f_loc, kind(f))
nf = int(nf_loc, kind(nf))
info = int(info_loc, kind(info))

contains

!--------------------------------------------------------------------------------------------------!
! This subroutine defines `calfun` using the C function pointer with an internal subroutine.
! This allows to avoid passing the C function pointer by a module variable, which is thread-unsafe.
! A possible security downside is that the compiler must allow for an executable stack.
!--------------------------------------------------------------------------------------------------!
subroutine calfun(x_sub, f_sub)
use, non_intrinsic :: consts_mod, only : RP
use, intrinsic :: iso_c_binding, only : C_DOUBLE
implicit none
real(RP), intent(in) :: x_sub(:)
real(RP), intent(out) :: f_sub

! Local variables
real(C_DOUBLE) :: x_sub_loc(size(x_sub))
real(C_DOUBLE) :: f_sub_loc

! Read the inputs and convert them to the types specified in COBJ
x_sub_loc = real(x_sub, kind(x_sub_loc))

! Call the C objective function
call obj_ptr(x_sub_loc, f_sub_loc, data_ptr)

! Write the output
f_sub = real(f_sub_loc, kind(f_sub))

end subroutine calfun

! We name some variables _sub to avoid masking the parent variables
subroutine callback_fcn(x_sub, f_sub, nf_sub, tr, cstrv, nlconstr, terminate)
use, non_intrinsic :: consts_mod, only : RP, IK
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_INT, C_BOOL
implicit none
real(RP), intent(in) :: x_sub(:)
real(RP), intent(in) :: f_sub
integer(IK), intent(in) :: nf_sub
integer(IK), intent(in) :: tr
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: nlconstr(:)
logical, intent(out) :: terminate

! Local variables
integer(C_INT) :: n_sub_loc
real(C_DOUBLE) :: x_sub_loc(size(x_sub))
real(C_DOUBLE) :: f_sub_loc
integer(C_INT) :: nf_sub_loc
integer(C_INT) :: tr_loc
real(C_DOUBLE) :: cstrv_loc
integer(C_INT) :: m_nlconstr
real(C_DOUBLE) :: nlconstr_loc(size(nlconstr))
logical(C_BOOL) :: terminate_loc

! Read the inputs and convert them to the types specified in CCALLBACK
n_sub_loc = size(x_sub)
x_sub_loc = real(x_sub, kind(x_sub_loc))
f_sub_loc = real(f_sub, kind(f_sub_loc))
nf_sub_loc = int(nf_sub, kind(nf_sub_loc))
tr_loc = int(tr, kind(tr_loc))
cstrv_loc = real(cstrv, kind(cstrv_loc))
m_nlconstr = size(nlconstr)
nlconstr_loc = real(nlconstr, kind(nlconstr_loc))

! Call the C objective function
call cb_ptr(n_sub_loc, x_sub_loc, f_sub_loc, nf_sub_loc, tr_loc, cstrv_loc, m_nlconstr, nlconstr_loc, terminate_loc)

! Write the output
terminate = logical(terminate_loc, kind(terminate))

end subroutine callback_fcn


end subroutine uobyqa_c


end module uobyqa_c_mod
