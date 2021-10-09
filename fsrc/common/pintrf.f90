module pintrf_mod
!--------------------------------------------------------------------------------------------------!
! This is a module specifying the abstract interfaces FUN and FUNCON. FUN evaluates the objective
! function for unconstrained, bound constrained, and linearly constrained problems; FUNCON evaluates
! the objective and constraint functions for nonlinearly constrained problems.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Saturday, October 09, 2021 PM11:01:46
!--------------------------------------------------------------------------------------------------!

!!!!!! Users must provide the implementation of FUN or FUNCON. !!!!!!

implicit none
private
public :: FUN, FUNCON

abstract interface
    subroutine FUN(x, f)
    use consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: f
    end subroutine FUN
end interface

abstract interface
    subroutine FUNCON(x, f, constr)
    use consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: f
    real(RP), intent(out) :: constr(:)
    end subroutine FUNCON
end interface

end module pintrf_mod
