! PROB_MOD is a module specifying the abstract interfaces FUNEVAL and 
! FCEVAL. FUNEVAL evaluates the objective function for unconstrained, 
! bound constrained, and linearly constrained problems; FCEVAL evaluates
! the objective function and constraint for nonlinearly constrained prolems.
!
! Coded by Zaikun Zhang in July 2020.
!
! Users must provide the implementation of CALFUN or CALCFC.

module prob_mod

implicit none
private
public :: FUNEVAL, FCEVAL 

abstract interface 
    subroutine FUNEVAL(x, f)
    use consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: f
    end subroutine FUNEVAL
end interface


abstract interface
    subroutine FCEVAL(x, f, con)
    use consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: f
    real(RP), intent(out) :: con(:)
    end subroutine FCEVAL 
end interface

end module prob_mod

