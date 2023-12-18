module cintrf_mod
!--------------------------------------------------------------------------------------------------!
! cintrf_mod provides interoperability with C objective functions
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: COBJ, COBJCON, CCALLBACK


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


end module cintrf_mod
