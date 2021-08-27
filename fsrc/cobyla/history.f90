module history_mod

implicit none
private
public :: savehist, savefilt

contains


subroutine savehist(nf, con, cstrv, f, x, chist, conhist, fhist, xhist)
! This subroutine saves X, F, CON, and CSTRV into XHIST, FHIST, CONHIST, and CHIST respectively.
use consts_mod, only : RP, IK, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize
implicit none

! Inputs
integer(IK), intent(in) :: nf
real(RP), intent(in) :: con(:)
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
character(len=*), parameter :: srname = 'SAVEHIST'

! Get and verify the sizes.
m = size(con)
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
! CON, and CSTRV at the position indexed by MOD(NF - 1, MAXHIST) + 1. Before the solver terminates,
! the history will be re-ordered so that the information is in the chronological order.
if (maxxhist > 0) then
    xhist(:, mod(nf - 1_IK, maxhist) + 1_IK) = x
end if
if (maxfhist > 0) then
    fhist(mod(nf - 1_IK, maxhist) + 1_IK) = f
end if
if (maxconhist > 0) then
    conhist(:, mod(nf - 1_IK, maxhist) + 1_IK) = con
end if
if (maxchist > 0) then
    chist(mod(nf - 1_IK, maxhist) + 1_IK) = cstrv
end if

end subroutine savehist


subroutine savefilt(con, cstrv, ctol, f, x, nfilt, cfilt, confilt, ffilt, xfilt)
! This subroutine saves XDROP in XSAV and DATDROP in DATSAV, unless a vector in XSAV(:, 1:NFILT) is
! better than XDROP. If XDROP is better than some vectors in XSAV(:, 1:NFILT), then these vectors
! will be removed. If XDROP is not better than any of XSAV(:, 1:NFILT) but NFILT=NFILTMAX, then we
! remove XSAV(:,1), which is the oldest vector in XSAV(:, 1:NFILT).
!
! When COBYLA calls this subroutine, XDROP is a vector to be "dropped", and  DATDROP contains its
! function/constraint information (constraint value in the first M entries, DATDROP(M+1) = F(XDROP),
! and DATDROP(M+2) = CSTRV(X)). XSAV and DATSAV save at most NFILTMAX vectors "dropped" by COBYLB
! and their function/constraint information. Only XSAV(:, 1:NFILT) and DATSAV(:, 1:NFILT) contains
! such vectors, while XSAV(:, NFILT+1:NFILTMAX) and DATSAV(:, NFILT+1:NFILTMAX) are not initialized yet.
!
! Note: We decide whether X is better than the function/constraint of Y by the ISBETTER function.

use consts_mod, only : RP, IK, DEBUGGING, SRNLEN
use debug_mod, only : errstop, verisize
use memory_mod, only : safealloc
use selectx_mod, only : isbetter

implicit none

! Inputs
real(RP), intent(in) :: con(:)  ! M
real(RP), intent(in) :: cstrv
real(RP), intent(in) :: ctol
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)  ! N

! In-outputs
integer(IK), intent(inout) :: nfilt
real(RP), intent(inout) :: cfilt(:)  ! NFILTMAX
real(RP), intent(inout) :: confilt(:, :)  ! (M, NFILTMAX)
real(RP), intent(inout) :: ffilt(:)  ! NFILTMAX
real(RP), intent(inout) :: xfilt(:, :) ! (N, NFILTMAX)

! Local variables
integer(IK) :: m
integer(IK) :: n
integer(IK) :: nfiltmax
integer(IK) :: i
integer(IK), allocatable :: index_to_keep(:)
logical :: better(nfilt)
logical :: keep(nfilt)
character(len=*), parameter :: srname = 'SAVEFILT'

m = size(con)
n = size(x)
nfiltmax = size(ffilt)

if (DEBUGGING) then
    if (n == 0) then
        call errstop(srname, 'SIZE(X) == 0')
    end if
    if (nfiltmax == 0) then
        call errstop(srname, 'SIZE(FFILT) == 0')
    end if
    call verisize(xfilt, n, nfiltmax)
    call verisize(confilt, m, nfiltmax)
    call verisize(cfilt, nfiltmax)
end if

! Return immediately if any column of XSAV is better than XDROP.
! BETTER is defined by the array constructor with an implicit do loop.
better = [(isbetter([ffilt(i), cfilt(i)], [f, cstrv], ctol), i=1, nfilt)]
if (any(better)) then
    return
end if

! Decide which columns of XSAV to keep. We use again the array constructor with an implicit do loop.
keep = [(.not. isbetter([f, cstrv], [ffilt(i), cfilt(i)], ctol), i=1, nfilt)]
! If X is not better than any column of XSAV, then we remove the first (oldest) column of XSAV.
if (count(keep) == nfiltmax) then
    keep(1) = .false.
end if

!--------------------------------------------------!
!----The SAFEALLOC line is removable in F2003.-----!
call safealloc(index_to_keep, nfilt)
!--------------------------------------------------!
index_to_keep = pack([(i, i=1, nfilt)], mask=keep)
nfilt = count(keep) + 1
xfilt(:, 1:nfilt - 1) = xfilt(:, index_to_keep)
ffilt(1:nfilt - 1) = ffilt(index_to_keep)
confilt(:, 1:nfilt - 1) = confilt(:, index_to_keep)
cfilt(1:nfilt - 1) = cfilt(index_to_keep)
! F2003 automatically deallocate local ALLOCATABLE variables at exit, yet
! we prefer to deallocate them immediately when they finish their jobs.
deallocate (index_to_keep)
xfilt(:, nfilt) = x
ffilt(nfilt) = f
confilt(:, nfilt) = con
cfilt(nfilt) = cstrv

end subroutine savefilt


end module history_mod
