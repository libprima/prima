module history_mod

implicit none
private
public :: savehist

interface savehist
    module procedure savehist_unc, savehist_nlc
end interface savehist


contains


subroutine savehist_unc(nf, f, x, fhist, xhist)
! This subroutine saves X and F into XHIST and FHIST respectively.
use consts_mod, only : RP, IK, DEBUGGING
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
character(len=*), parameter :: srname = 'SAVEHIST_UNC'

! Get and verify the sizes.
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))
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

end subroutine savehist_unc


subroutine savehist_nlc(nf, constr, cstrv, f, x, chist, conhist, fhist, xhist)
! This subroutine saves X, F, CONSTR, and CSTRV into XHIST, FHIST, CONHIST, and CHIST respectively.
use consts_mod, only : RP, IK, DEBUGGING
use debug_mod, only : errstop, verisize
implicit none

! Inputs
integer(IK), intent(in) :: nf
real(RP), intent(in) :: constr(:)
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)

! In-outputs
real(RP), intent(inout) :: chist(:)
real(RP), intent(inout) :: conhist(:, :)
real(RP), intent(inout) :: fhist(:)
real(RP), intent(inout) :: xhist(:, :)

! Local variables
integer(IK) :: m
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
character(len=*), parameter :: srname = 'SAVEHIST_NLC'

! Get and verify the sizes.
m = size(constr)
n = size(x)
maxxhist = size(xhist, 2)
maxfhist = size(fhist)
maxconhist = size(conhist, 2)
maxchist = size(chist)
maxhist = max(maxxhist, maxfhist, maxconhist, maxchist)
if (DEBUGGING) then
    if (n < 1) then
        call errstop(srname, 'SIZE(X) is invalid')
    end if
    if (maxxhist > 0) then
        call verisize(xhist, n, maxhist)
    end if
    if (maxfhist > 0) then
        call verisize(fhist, maxhist)
    end if
    if (maxconhist > 0) then
        call verisize(conhist, m, maxhist)
    end if
    if (maxchist > 0) then
        call verisize(chist, maxhist)
    end if
end if

! Save the history. Note that NF may exceed the maximal amount of history to save. We save X, F,
! CONSTR, and CSTRV at the position indexed by MOD(NF - 1, MAXHIST) + 1. Before the solver terminates,
! the history will be re-ordered so that the information is in the chronological order.
if (maxxhist > 0) then
    xhist(:, mod(nf - 1_IK, maxhist) + 1_IK) = x
end if
if (maxfhist > 0) then
    fhist(mod(nf - 1_IK, maxhist) + 1_IK) = f
end if
if (maxconhist > 0) then
    conhist(:, mod(nf - 1_IK, maxhist) + 1_IK) = constr
end if
if (maxchist > 0) then
    chist(mod(nf - 1_IK, maxhist) + 1_IK) = cstrv
end if

end subroutine savehist_nlc


end module history_mod
