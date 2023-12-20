module cobyla_c_mod
!--------------------------------------------------------------------------------------------------!
! cobyla_c_mod provides cobyla_c, a simplified interface to cobyla for interoperability with C
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!--------------------------------------------------------------------------------------------------!
implicit none
private
public :: cobyla_c


contains


subroutine cobyla_c(m_nlcon, cobjcon_ptr, data_ptr, n, x, f, cstrv, nlconstr, m_ineq, Aineq, bineq, m_eq, Aeq, beq, &
    & xl, xu, f0, nlconstr0, nf, rhobeg, rhoend, ftarget, maxfun, iprint, callback_ptr, info) bind(C)
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_INT, C_FUNPTR, C_PTR, C_ASSOCIATED, C_F_PROCPOINTER
use, non_intrinsic :: cintrf_mod, only : COBJCON, CCALLBACK
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: cobyla_mod, only : cobyla
implicit none

! Compulsory arguments
integer(C_INT), intent(in), value :: m_nlcon
type(C_FUNPTR), intent(in), value :: cobjcon_ptr
type(C_PTR), intent(in), value :: data_ptr
integer(C_INT), intent(in), value :: n
! We cannot use assumed-shape arrays for C interoperability
real(C_DOUBLE), intent(inout) :: x(n)
real(C_DOUBLE), intent(out) :: f
real(C_DOUBLE), intent(out) :: cstrv
real(C_DOUBLE), intent(out) :: nlconstr(m_nlcon)
integer(C_INT), intent(in), value :: m_ineq
real(C_DOUBLE), intent(in) :: Aineq(n, m_ineq)
real(C_DOUBLE), intent(in) :: bineq(m_ineq)
integer(C_INT), intent(in), value :: m_eq
real(C_DOUBLE), intent(in) :: Aeq(n, m_eq)
real(C_DOUBLE), intent(in) :: beq(m_eq)
real(C_DOUBLE), intent(in) :: xl(n)
real(C_DOUBLE), intent(in) :: xu(n)
real(C_DOUBLE), intent(in), value :: f0
real(C_DOUBLE), intent(in) :: nlconstr0(m_nlcon)
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
integer(IK) :: m_nlcon_loc
integer(IK) :: maxfun_loc
integer(IK) :: nf_loc
real(RP) :: Aineq_loc(m_ineq, n)
real(RP) :: bineq_loc(m_ineq)
real(RP) :: Aeq_loc(m_eq, n)
real(RP) :: beq_loc(m_eq)
real(RP) :: cstrv_loc
real(RP) :: nlconstr_loc(m_nlcon)
real(RP) :: f_loc
real(RP) :: rhobeg_loc
real(RP) :: rhoend_loc
real(RP) :: ftarget_loc
real(RP) :: x_loc(n)
real(RP) :: xl_loc(n)
real(RP) :: xu_loc(n)
real(RP) :: f0_loc
real(RP) :: nlconstr0_loc(m_nlcon)
! The initialization to null is necessary to avoid a bug with the newer Intel compiler ifx.
! See details here: https://fortran-lang.discourse.group/t/strange-issue-with-ifx-compiler-and-assume-recursion/7013
! The bug was observed in all versions of ifx up to 2024.0.1. Once this bug is fixed we should remove the
! initialization to null because it implies the 'save' attribute, which is undesirable.
procedure(COBJCON), pointer :: objcon_ptr => null()
procedure(CCALLBACK), pointer :: cb_ptr => null()

! Read the inputs and convert them to the Fortran side types
! Note that `transpose` is needed when reading 2D arrays, since they are stored in the row-major
! order in c but column-major in Fortran.
x_loc = real(x, kind(x_loc))
Aineq_loc = real(transpose(Aineq), kind(Aineq_loc))
bineq_loc = real(bineq, kind(bineq_loc))
Aeq_loc = real(transpose(Aeq), kind(Aeq_loc))
beq_loc = real(beq, kind(beq_loc))
xl_loc = real(xl, kind(xl_loc))
xu_loc = real(xu, kind(xu_loc))
f0_loc = real(f0, kind(f0_loc))
nlconstr0_loc = real(nlconstr0, kind(nlconstr0_loc))
rhobeg_loc = real(rhobeg, kind(rhobeg_loc))
rhoend_loc = real(rhoend, kind(rhoend_loc))
ftarget_loc = real(ftarget, kind(ftarget_loc))
maxfun_loc = int(maxfun, kind(maxfun_loc))
iprint_loc = int(iprint, kind(iprint_loc))
m_nlcon_loc = int(m_nlcon, kind(m_nlcon_loc))
call C_F_PROCPOINTER(cobjcon_ptr, objcon_ptr)

! Call the Fortran code
if (C_ASSOCIATED(callback_ptr)) then
    ! If a C callback function is provided, we convert it to a Fortran procedure pointer and capture
    ! that pointer in the closure below.
    call C_F_PROCPOINTER(callback_ptr, cb_ptr)
    ! We then provide the closure to the algorithm.
    call cobyla(calcfc, m_nlcon_loc, x_loc, f_loc, cstrv=cstrv_loc, nlconstr=nlconstr_loc, Aineq=Aineq_loc, &
        & bineq=bineq_loc, Aeq=Aeq_loc, beq=beq_loc, xl=xl_loc, xu=xu_loc, f0=f0_loc, nlconstr0=nlconstr0_loc, &
        & nf=nf_loc, rhobeg=rhobeg_loc, rhoend=rhoend_loc, ftarget=ftarget_loc, maxfun=maxfun_loc, &
        & iprint=iprint_loc, callback_fcn=callback_fcn, info=info_loc)
else
    call cobyla(calcfc, m_nlcon_loc, x_loc, f_loc, cstrv=cstrv_loc, nlconstr=nlconstr_loc, Aineq=Aineq_loc, &
        & bineq=bineq_loc, Aeq=Aeq_loc, beq=beq_loc, xl=xl_loc, xu=xu_loc, f0=f0_loc, nlconstr0=nlconstr0_loc, &
        & nf=nf_loc, rhobeg=rhobeg_loc, rhoend=rhoend_loc, ftarget=ftarget_loc, maxfun=maxfun_loc, &
        & iprint=iprint_loc, info=info_loc)
end if

! Write the outputs
x = real(x_loc, kind(x))
f = real(f_loc, kind(f))
cstrv = real(cstrv_loc, kind(cstrv))
nf = int(nf_loc, kind(nf))
info = int(info_loc, kind(info))
nlconstr = real(nlconstr_loc, kind(nlconstr))

contains

!--------------------------------------------------------------------------------------------------!
! This subroutine defines `calcfc` using the C function pointer with an internal subroutine.
! This allows to avoid passing the C function pointer by a module variable, which is thread-unsafe.
! A possible security downside is that the compiler must allow for an executable stack.
! This subroutine is identical across 4 out of 5 algorithms; COBYLA requires a slightly different
! signature.
!--------------------------------------------------------------------------------------------------!
subroutine calcfc(x_sub, f_sub, constr_sub)
use, intrinsic :: iso_c_binding, only : C_DOUBLE
use, non_intrinsic :: consts_mod, only : RP
implicit none
real(RP), intent(in) :: x_sub(:) ! We name some variables _sub to avoid masking the parent variables
real(RP), intent(out) :: f_sub
real(RP), intent(out) :: constr_sub(:)

! Local variables
real(C_DOUBLE) :: x_sub_loc(size(x_sub))
real(C_DOUBLE) :: f_sub_loc
real(C_DOUBLE) :: constr_sub_loc(size(constr_sub))

! Read the inputs and convert them to the types specified in COBJCON
x_sub_loc = real(x_sub, kind(x_sub_loc))

! Call the C objective function
call objcon_ptr(x_sub_loc, f_sub_loc, constr_sub_loc, data_ptr)

! Write the output
f_sub = real(f_sub_loc, kind(f_sub))
constr_sub = real(constr_sub_loc, kind(constr_sub))

end subroutine calcfc


!--------------------------------------------------------------------------------------------------!
! This subroutine defines `callback_fcn` using the C function pointer with an internal subroutine.
! This allows to avoid passing the C function pointer by a module variable, which is thread-unsafe.
! A possible security downside is that the compiler must allow for an executable stack.
! This subroutine is identical across all 5 algorithms.
!--------------------------------------------------------------------------------------------------!
subroutine callback_fcn(x_sub, f_sub, nf_sub, tr, cstrv_sub, nlconstr_sub, terminate)
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_INT, C_BOOL
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

! Call the C callback function
call cb_ptr(n_sub_loc, x_sub_loc, f_sub_loc, nf_sub_loc, tr_loc, cstrv_sub_loc, m_nlconstr, nlconstr_sub_loc, terminate_loc)

! Write the output
if ( present(terminate) ) then
    terminate = logical(terminate_loc, kind(terminate))
end if

! Deallocate resources
if (allocated(nlconstr_sub_loc)) deallocate(nlconstr_sub_loc)

end subroutine callback_fcn


end subroutine cobyla_c


end module cobyla_c_mod
