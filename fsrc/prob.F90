! PROB_MOD is a module specifying the interfaces for CALFUN and CALCFC.
! CALFUN evaluates the objective function for unconstrained, bound
! constrained, and linearly constrained problems; CALCFC evaluates the
! objective function and constraint for nonlinearly constrained prolems.
!
! Coded by Zaikun Zhang in July 2020.
!
! Users must provide the implementation of CALFUN or CALCFC.

module prob_mod

implicit none
private

public :: calfun, calcfc

interface 
    subroutine calfun(x, f)
    use consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: f
    end subroutine calfun
end interface


interface
    subroutine calcfc(x, f, con)
    use consts_mod, only : RP
    implicit none
    real(RP), intent(in) :: x(:)
    real(RP), intent(out) :: f
    real(RP), intent(out) :: con(:)
    end subroutine calcfc
end interface

end module prob_mod
