module pintrf_mod
!--------------------------------------------------------------------------------------------------!
! This is a module specifying the abstract interfaces OBJ, OBJCON, and CALLBACK. OBJ evaluates the
! objective function for unconstrained, bound constrained, and linearly constrained problems; OBJCON
! evaluates the objective and constraint functions for nonlinearly constrained problems; CALLBACK
! is a callback function that is called after each iteration of the solvers to report the progress
! and optionally request termination.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Friday, December 22, 2023 PM01:23:44
!--------------------------------------------------------------------------------------------------!

!!!!!! Users must provide the implementation of OBJ or OBJCON. !!!!!!

implicit none
private
public :: OBJ, OBJCON, CALLBACK


abstract interface

    subroutine OBJ(x, f)
    use consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: f
    end subroutine OBJ


    subroutine OBJCON(x, f, constr)
    use consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: f
    real(RP), intent(out) :: constr(:)
    end subroutine OBJCON


    subroutine CALLBACK(x, f, nf, tr, cstrv, nlconstr, terminate)
    use consts_mod, only : RP, IK
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(in) :: f
    integer(IK), intent(in) :: nf
    integer(IK), intent(in) :: tr
    real(RP), intent(in), optional :: cstrv
    real(RP), intent(in), optional :: nlconstr(:)
    logical, intent(out), optional :: terminate
    end subroutine CALLBACK

end interface

end module pintrf_mod
