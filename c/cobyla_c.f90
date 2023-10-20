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


subroutine cobyla_c(m_nlcon, cobjcon_ptr, data_ptr, n, x, f, cstrv, nlconstr, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu, &
    & nf, rhobeg, rhoend, ftarget, maxfun, iprint, info) bind(C)
use, intrinsic :: iso_c_binding, only : C_DOUBLE, C_INT, C_FUNPTR, C_PTR
use, non_intrinsic :: cintrf_mod, only : COBJCON
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
integer(C_INT), intent(out) :: nf
real(C_DOUBLE), intent(in), value :: rhobeg
real(C_DOUBLE), intent(in), value :: rhoend
real(C_DOUBLE), intent(in), value :: ftarget
integer(C_INT), intent(in), value :: maxfun
integer(C_INT), intent(in), value :: iprint
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
rhobeg_loc = real(rhobeg, kind(rhobeg_loc))
rhoend_loc = real(rhoend, kind(rhoend_loc))
ftarget_loc = real(ftarget, kind(ftarget_loc))
maxfun_loc = int(maxfun, kind(maxfun_loc))
iprint_loc = int(iprint, kind(iprint_loc))
m_nlcon_loc = int(m_nlcon, kind(m_nlcon_loc))

! Call the Fortran code
call cobyla(calcfc, m_nlcon_loc, x_loc, f_loc, cstrv=cstrv_loc, nlconstr=nlconstr_loc, &
    & Aineq=Aineq_loc, bineq=bineq_loc, Aeq=Aeq_loc, beq=beq_loc, &
    & xl=xl_loc, xu=xu_loc, nf=nf_loc, &
    & rhobeg=rhobeg_loc, rhoend=rhoend_loc, ftarget=ftarget_loc, maxfun=maxfun_loc, &
    & iprint=iprint_loc, info=info_loc)

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
!--------------------------------------------------------------------------------------------------!
subroutine calcfc(x_sub, f_sub, constr_sub)
use, non_intrinsic :: consts_mod, only : RP
use, non_intrinsic :: cintrf_mod, only : evalcobjcon
implicit none
real(RP), intent(in) :: x_sub(:)
real(RP), intent(out) :: f_sub
real(RP), intent(out) :: constr_sub(:)
call evalcobjcon(cobjcon_ptr, data_ptr, x_sub, f_sub, constr_sub)
end subroutine calcfc

end subroutine cobyla_c


end module cobyla_c_mod
