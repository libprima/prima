module newuoa_c_mod
!--------------------------------------------------------------------------------------------------!
! newuoa_c_mod provides newuoa_c, a simplified interface to cobyla for interoperability with C
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!--------------------------------------------------------------------------------------------------!
implicit none
private
public :: newuoa_c


contains


subroutine newuoa_c(cobj_ptr, data_ptr, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, callback_ptr, info) bind(C)
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_INT, C_FUNPTR, C_PTR, C_ASSOCIATED
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: newuoa_mod, only : newuoa
implicit none

! Compulsory arguments
type(C_FUNPTR), intent(IN), value :: cobj_ptr
type(C_PTR), intent(in) :: data_ptr
integer(C_INT), intent(in), value :: n
! We cannot use assumed-shape arrays for C interoperability
real(C_DOUBLE), intent(inout) :: x(n)
real(C_DOUBLE), intent(out) :: f
integer(C_INT), intent(out) :: nf
real(C_DOUBLE), intent(in), value :: rhobeg
real(C_DOUBLE), intent(in), value :: rhoend
real(C_DOUBLE), intent(in), value :: ftarget
integer(C_INT), intent(in), value :: maxfun
integer(C_INT), intent(in), value :: npt
integer(C_INT), intent(in), value :: iprint
type(C_FUNPTR), intent(in), value :: callback_ptr
integer(C_INT), intent(out) :: info

! Local variables
integer(IK) :: info_loc
integer(IK) :: iprint_loc
integer(IK) :: maxfun_loc
integer(IK) :: npt_loc
integer(IK) :: nf_loc
real(RP) :: f_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP) :: ftarget_loc
real(RP) :: x_loc(n)

! Read the inputs and convert them to the Fortran side types
x_loc = real(x, kind(x_loc))
rhobeg_loc = real(rhobeg, kind(rhobeg_loc))
rhoend_loc = real(rhoend, kind(rhoend_loc))
ftarget_loc = real(ftarget, kind(ftarget_loc))
maxfun_loc = int(maxfun, kind(maxfun_loc))
npt_loc = int(npt, kind(npt_loc))
iprint_loc = int(iprint, kind(iprint_loc))

! Call the Fortran code
if (C_ASSOCIATED(callback_ptr)) then
    ! If a C callback function is provided, we capture the callback_ptr for use in the closure below,
    ! and then we pass the closure to the Fortran code.
    call newuoa(calfun, x_loc, f_loc, nf=nf_loc, rhobeg=rhobeg_loc, rhoend=rhoend_loc, ftarget=ftarget_loc, &
        & maxfun=maxfun_loc, npt=npt_loc, iprint=iprint_loc, callback_fcn=callback_fcn, info=info_loc)
else
    call newuoa(calfun, x_loc, f_loc, nf=nf_loc, rhobeg=rhobeg_loc, rhoend=rhoend_loc, ftarget=ftarget_loc, &
        & maxfun=maxfun_loc, npt=npt_loc, iprint=iprint_loc, info=info_loc)
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
! This subroutine is identical across 4/5 algorithms; COBYLA requires a slightly different signature.
!--------------------------------------------------------------------------------------------------!
subroutine calfun(x_sub, f_sub)
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_F_PROCPOINTER
use, non_intrinsic :: cintrf_mod, only : COBJ
use, non_intrinsic :: consts_mod, only : RP
implicit none
real(RP), intent(in) :: x_sub(:) ! We name some variables _sub to avoid masking the parent variables
real(RP), intent(out) :: f_sub

! Local variables
real(C_DOUBLE) :: x_sub_loc(size(x_sub))
real(C_DOUBLE) :: f_sub_loc
procedure(COBJ), pointer :: obj_ptr

! Read the inputs and convert them to the types specified in COBJ
x_sub_loc = real(x_sub, kind(x_sub_loc))

! The intel compiler insists that we convert the C function pointer to a Fortran function pointer
! here as opposed to within the parent function, otherwise it segfaults.
call C_F_PROCPOINTER(cobj_ptr, obj_ptr)

! Call the C objective function
call obj_ptr(x_sub_loc, f_sub_loc, data_ptr)

! Write the output
f_sub = real(f_sub_loc, kind(f_sub))

end subroutine calfun


!--------------------------------------------------------------------------------------------------!
! This subroutine defines `callback_fcn` using the C function pointer with an internal subroutine.
! This allows to avoid passing the C function pointer by a module variable, which is thread-unsafe.
! A possible security downside is that the compiler must allow for an executable stack.
! This subroutine is identical across all 5 algorithms.
!--------------------------------------------------------------------------------------------------!
subroutine callback_fcn(x_sub, f_sub, nf_sub, tr, cstrv_sub, nlconstr_sub, terminate)
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_INT, C_BOOL, C_F_PROCPOINTER
use, non_intrinsic :: cintrf_mod, only : CCALLBACK
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: memory_mod, only : safealloc
implicit none
real(RP), intent(in) :: x_sub(:) ! We name some variables _sub to avoid masking the parent variables
real(RP), intent(in) :: f_sub
integer(IK), intent(in) :: nf_sub
integer(IK), intent(in) :: tr
real(RP), intent(in), optional :: cstrv_sub
real(RP), intent(in), optional :: nlconstr_sub(:)
logical, intent(out), optional :: terminate

! Local variables
integer(C_INT) :: n_sub_loc
real(C_DOUBLE) :: x_sub_loc(size(x_sub))
real(C_DOUBLE) :: f_sub_loc
integer(C_INT) :: nf_sub_loc
integer(C_INT) :: tr_loc
real(C_DOUBLE) :: cstrv_sub_loc
integer(C_INT) :: m_nlconstr
real(C_DOUBLE), allocatable :: nlconstr_sub_loc(:)
logical(C_BOOL) :: terminate_loc
procedure(CCALLBACK), pointer :: cb_ptr

! Read the inputs and convert them to the types specified in CCALLBACK
n_sub_loc = size(x_sub)
x_sub_loc = real(x_sub, kind(x_sub_loc))
f_sub_loc = real(f_sub, kind(f_sub_loc))
nf_sub_loc = int(nf_sub, kind(nf_sub_loc))
tr_loc = int(tr, kind(tr_loc))

! Set the constraint violation to a sensible default value if it is not provided.
if (present(cstrv_sub)) then
    cstrv_sub_loc = real(cstrv_sub, kind(cstrv_sub_loc))
else
    cstrv_sub_loc = 0.0_C_DOUBLE
end if

! Set the nonlinear constraints to a sensible default value if it is not provided.
if (present(nlconstr_sub)) then
    m_nlconstr = int(size(nlconstr_sub), C_INT)
    call safealloc(nlconstr_sub_loc, int(m_nlconstr, IK))
    nlconstr_sub_loc = real(nlconstr_sub, kind(nlconstr_sub_loc))
else
    m_nlconstr = 0_C_INT
    nlconstr_sub_loc = [real(C_DOUBLE) ::]
end if

! As above, the intel compiler insists on doing this conversion here, every time, as opposed to
! within the parent function, once.
call C_F_PROCPOINTER(callback_ptr, cb_ptr)

! Call the C objective function
call cb_ptr(n_sub_loc, x_sub_loc, f_sub_loc, nf_sub_loc, tr_loc, cstrv_sub_loc, m_nlconstr, nlconstr_sub_loc, terminate_loc)

! Write the output
if ( present(terminate) ) then
    terminate = logical(terminate_loc, kind(terminate))
end if

! Deallocate resources
if (allocated(nlconstr_sub_loc)) deallocate(nlconstr_sub_loc)

end subroutine callback_fcn

end subroutine newuoa_c


end module newuoa_c_mod
