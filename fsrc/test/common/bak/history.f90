module history_mod

implicit none
private
public :: savehist, rangehist

interface savehist
    module procedure savehist_unc, savehist_nlc
end interface savehist

interface rangehist
    module procedure rangehist_unc, rangehist_nlc
end interface rangehist


contains


subroutine savehist_unc(nf, f, x, fhist, xhist)
! This subroutine saves X and F into XHIST and FHIST respectively.
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : errstop, verisize
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
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : errstop, verisize
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
m = int(size(constr), kind(m))
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxconhist = int(size(conhist, 2), kind(maxconhist))
maxchist = int(size(chist), kind(maxchist))
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


subroutine rangehist_unc(nf, fhist, xhist)
! This subroutine arranges FHIST and XHIST in the chronological order.
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
integer(IK), intent(in) :: nf

! In-outputs
real(RP), intent(inout) :: fhist(:)
real(RP), intent(inout) :: xhist(:, :)

! Local variables
integer(IK) :: khist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
character(len=*), parameter :: srname = 'RANGEHIST_UNC'

! Get and verify the sizes.
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))
if (DEBUGGING) then
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(xhist, 1) >= 1 .and. maxxhist * (maxxhist - maxhist) == 0, &
         & 'SIZE(XHIST, 1) >= 1, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
end if

if (maxxhist > 0 .and. maxxhist < nf) then
    khist = mod(nf - 1_IK, maxxhist) + 1_IK
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
    ! N.B.:
    ! 1. The result of the array constructor is always a rank-1 array (e.g., vector), no matter what
    ! elements are used to construct the array.
    ! 2. The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
    ! order of Fortran arrays.
    ! 3. In MATLAB, `xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)]` does the same thing.
end if
if (maxfhist > 0 .and. maxfhist < nf) then
    khist = mod(nf - 1_IK, maxfhist) + 1_IK
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
end if

end subroutine rangehist_unc


subroutine rangehist_nlc(nf, chist, conhist, fhist, xhist)
! This subroutine arranges FHIST and XHIST, CONHIST, and CHIST in the chronological order.
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
integer(IK), intent(in) :: nf

! In-outputs
real(RP), intent(inout) :: chist(:)
real(RP), intent(inout) :: conhist(:, :)
real(RP), intent(inout) :: fhist(:)
real(RP), intent(inout) :: xhist(:, :)

! Local variables
integer(IK) :: khist
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
character(len=*), parameter :: srname = 'RANGEHIST_NLC'

! Get and verify the sizes.
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxconhist = int(size(conhist, 2), kind(maxconhist))
maxchist = int(size(chist), kind(maxchist))
maxhist = max(maxxhist, maxfhist, maxconhist, maxchist)
if (DEBUGGING) then
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
    call assert(maxconhist * (maxconhist - maxhist) == 0, &
         & 'SIZE(CONHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(xhist, 1) >= 1 .and. maxxhist * (maxxhist - maxhist) == 0, &
         & 'SIZE(XHIST, 1) >= 1, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
end if

if (maxxhist > 0 .and. maxxhist < nf) then
    khist = mod(nf - 1_IK, maxxhist) + 1_IK
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
    ! N.B.:
    ! 1. The result of the array constructor is always a rank-1 array (e.g., vector), no matter what
    ! elements are used to construct the array.
    ! 2. The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
    ! order of Fortran arrays.
    ! 3. In MATLAB, `xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)]` does the same thing.
end if
if (maxfhist > 0 .and. maxfhist < nf) then
    khist = mod(nf - 1_IK, maxfhist) + 1_IK
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
end if
if (maxconhist > 0 .and. maxconhist < nf) then
    khist = mod(nf - 1_IK, maxconhist) + 1_IK
    conhist = reshape([conhist(:, khist + 1:maxconhist), conhist(:, 1:khist)], shape(conhist))
end if
if (maxchist > 0 .and. maxchist < nf) then
    khist = mod(nf - 1_IK, maxchist) + 1_IK
    chist = [chist(khist + 1:maxchist), chist(1:khist)]
end if

end subroutine rangehist_nlc


end module history_mod
