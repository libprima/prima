module history_mod

implicit none
private
public :: savehist

contains


subroutine savehist(nf, f, x, fhist, xhist)
! This subroutine saves X and F into XHIST and FHIST respectively.
use consts_mod, only : RP, IK, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize
implicit none

! Inputs
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)

! In-outputs
real(RP), intent(inout) :: fhist(:)
real(RP), intent(inout) :: xhist(:, :)

! Local variables
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
character(len=SRNLEN), parameter :: srname = 'SAVEHIST'

! Get and verify the sizes.
n = size(x)
maxxhist = size(xhist, 2)
maxfhist = size(fhist)
maxhist = max(maxxhist, maxfhist)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(X) < 1')
    end if
    if (maxxhist > 0) then
        call verisize(xhist, n, maxhist)
    end if
    if (maxfhist > 0) then
        call verisize(fhist, maxhist)
    end if
end if

! Save the history. Note that NF may exceed the maximal amount of history to save. We save X and F
! at the position indexed by MOD(NF - 1, MAXHIST) + 1. Before the solver terminates, the history
! will be re-ordered so that the information is in the chronological order.
if (maxxhist > 0) then
    xhist(:, mod(nf - 1_IK, maxhist) + 1_IK) = x
end if
if (maxfhist > 0) then
    fhist(mod(nf - 1_IK, maxhist) + 1_IK) = f
end if

end subroutine savehist


end module history_mod
