module history_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines that handle the X/F/C histories of the solver, taking into
! account that MAXHIST may be smaller than NF.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Monday, March 14, 2022 PM02:23:43
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: prehist
public :: savehist
public :: rangehist


contains


subroutine prehist(maxhist, n, output_xhist, xhist, output_fhist, fhist, output_chist, chist, m, output_conhist, conhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine revises MAXHIST according to MAXMEMORY, and allocates memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! XHIST = NaN(N, MAXFUN), FHIST = NaN(1, MAXFUN), CHIST = NaN(1, MAXFUN), CONHIST = NaN(M, MAXFUN),
! if they are requested; replace MAXFUN with 0 for the history that is not requested.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, MAXMEMORY, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : int
use, non_intrinsic :: memory_mod, only : safealloc, cstyle_sizeof
implicit none

! Inputs
integer(IK), intent(in) :: n
integer(IK), intent(in), optional :: m
logical, intent(in) :: output_fhist
logical, intent(in) :: output_xhist
logical, intent(in), optional :: output_chist
logical, intent(in), optional :: output_conhist

! In-outputs
integer(IK), intent(inout) :: maxhist

! Outputs
real(RP), intent(out), allocatable :: fhist(:)
real(RP), intent(out), allocatable :: xhist(:, :)
real(RP), intent(out), optional, allocatable :: chist(:)
real(RP), intent(out), optional, allocatable :: conhist(:, :)

! Local variables
character(len=*), parameter :: srname = 'PREHIST'
integer(IK) :: maxhist_in
integer(IK) :: unit_memo

! Preconditions
if (DEBUGGING) then
    call assert(maxhist >= 0, 'MAXHIST >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
    if (present(m)) then
        call assert(m >= 0, 'M >= 0', srname)
    end if
    call assert(present(output_chist) .eqv. present(chist), &
        & 'OUTPUT_CHIST and CHIST are both present or both absent', srname)
    call assert((present(m) .eqv. present(conhist)) .and. (present(output_conhist) .eqv. present(conhist)), &
        & 'M, OUTPUT_CONHIST, and CONHIST are all present or all absent', srname)
end if

!====================!
! Calculation starts !
!====================!

! Save the input value of MAXHIST for debugging.
maxhist_in = maxhist

! Revise MAXHIST according to MAXMEMORY, i.e., the maximal memory allowed for the history.
unit_memo = int(output_xhist) * n + int(output_fhist)
if (present(output_chist) .and. present(chist)) then
    unit_memo = unit_memo + int(output_chist)
end if
if (present(m) .and. present(output_conhist) .and. present(conhist)) then
    unit_memo = unit_memo + int(output_conhist) * m
end if
unit_memo = unit_memo * cstyle_sizeof(0.0_RP)
if (unit_memo <= 0) then  ! No output of history is requested
    maxhist = 0_IK
elseif (maxhist > MAXMEMORY / unit_memo) then
    maxhist = int(MAXMEMORY / unit_memo, kind(maxhist))
    ! We cannot simply set MAXHIST = MIN(MAXHIST, MAXMEMORY/UNIT_MEMO), as they may not have
    ! the same kind, and compilers may complain. We may convert them, but overflow may occur.
end if

call safealloc(xhist, n, maxhist * int(output_xhist))
call safealloc(fhist, maxhist * int(output_fhist))
! Even if OUTPUT_CHIST is FALSE, CHIST still needs to be allocated.
if (present(output_chist) .and. present(chist)) then
    call safealloc(chist, maxhist * int(output_chist))
end if
! Even if OUTPUT_CONHIST is FALSE, CONHIST still needs to be allocated.
if (present(m) .and. present(output_conhist) .and. present(conhist)) then
    call safealloc(conhist, m, maxhist * int(output_conhist))
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(maxhist >= 0 .and. maxhist <= maxhist_in, '0 <= MAXHIST <= MAXHIST_IN', srname)
    call assert(int(maxhist, kind(MAXMEMORY)) * int(unit_memo, kind(MAXMEMORY)) <= MAXMEMORY, &
        & 'The history will not take more memory than MAXMEMORY', srname)
    call assert(allocated(xhist), 'XHIST is allocated', srname)
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxhist * int(output_xhist), &
        & 'if XHIST is requested, then SIZE(XHIST) == [N, MAXHIST]; otherwise, SIZE(XHIST) == [N, 0]', srname)
    call assert(allocated(fhist), 'FHIST is allocated', srname)
    call assert(size(fhist) == maxhist * int(output_fhist), &
        & 'if FHIST is requested, then SIZE(FHIST) == MAXHIST; otherwise, SIZE(FHIST) == 0', srname)
    if (present(output_chist) .and. present(chist)) then
        call assert(allocated(chist), 'CHIST is allocated', srname)
        call assert(size(chist) == maxhist * int(output_chist), &
            & 'if CHIST is requested, then SIZE(CHIST) == MAXHIST; otherwise, SIZE(CHIST) == 0', srname)
    end if
    if (present(m) .and. present(output_conhist) .and. present(conhist)) then
        call assert(allocated(conhist), 'CONHIST is allocated', srname)
        call assert(size(conhist, 1) == m .and. size(conhist, 2) == maxhist * int(output_conhist), &
            & 'if CONHIST is requested, then SIZE(CONHIST) == [M, MAXHIST]; otherwise, SIZE(CONHIST) == [M, 0]', srname)
    end if
end if
end subroutine prehist


subroutine savehist(nf, x, xhist, f, fhist, cstrv, chist, constr, conhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine saves X, F, CSTRV, and CONSTR into XHIST, FHIST, CHIST, and CONHIST respectively.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf, is_neginf
implicit none

! Inputs
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)
real(RP), intent(in), optional :: constr(:)
real(RP), intent(in), optional :: cstrv

! In-outputs
real(RP), intent(inout) :: fhist(:)
real(RP), intent(inout) :: xhist(:, :)
real(RP), intent(inout), optional :: chist(:)
real(RP), intent(inout), optional :: conhist(:, :)

! Local variables
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
character(len=*), parameter :: srname = 'SAVEHIST'

! Sizes
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
if (present(chist) .and. present(cstrv)) then
    maxchist = int(size(chist), kind(maxchist))
else
    maxchist = 0_IK
end if
if (present(conhist) .and. present(constr)) then
    maxconhist = int(size(conhist, 2), kind(maxconhist))
else
    maxconhist = 0_IK
end if
maxhist = max(maxxhist, maxfhist, maxchist, maxconhist)

! Preconditions
if (DEBUGGING) then  ! Called after each function evaluation when debugging; can be expensive.
    ! Check the presence of CSTRV, CHIST, CONSTR, CONHIST.
    call assert(present(cstrv) .eqv. present(chist), 'CSTRV and CHIST are both present or both absent', srname)
    call assert(present(constr) .eqv. present(conhist), 'CONSTR and CONHIST are both present or both absent', srname)
    ! Check the size of X.
    call assert(size(x) >= 1, 'SIZE(X) >= 1', srname)
    ! Check the sizes of XHIST, FHIST, CONHIST, CHIST.
    call assert(size(xhist, 1) == size(x) .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == SIZE(X), SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
    if (present(constr) .and. present(conhist)) then
        call assert(size(conhist, 1) == size(constr) .and. maxconhist * (maxconhist - maxhist) == 0, &
            & 'SIZE(CONHIST, 1) == SIZE(CONSTR), SIZE(CONHIST, 2) == 0 or MAXNHIST', srname)
    end if
    ! Check the values of XHIST, FHIST, CHIST, CONHIST, up to the (NF - 1)th position.
    ! As long as this subroutine is called, XHIST contains only finite values.
    call assert(all(is_finite(xhist(:, 1:min(nf - 1_IK, maxxhist)))), 'XHIST is finite', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf - 1_IK, maxfhist))) .or. &
        & is_posinf(fhist(1:min(nf - 1_IK, maxfhist)))), 'FHIST does not contain NaN/+Inf', srname)
    !----------------------------------------------------------------------------------------------!
    ! The following test is not applicable to LINCOA.
    !!if (present(chist)) then
    !!    call assert(.not. any(chist(1:min(nf - 1_IK, maxchist)) < 0 .or. &
    !!        & is_nan(chist(1:min(nf - 1_IK, maxchist))) .or. is_posinf(chist(1:min(nf - 1_IK, maxchist)))), &
    !!        & 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
    !!end if
    !----------------------------------------------------------------------------------------------!
    if (present(conhist)) then
        call assert(.not. any(is_nan(conhist(:, 1:min(nf - 1_IK, maxconhist))) .or. &
            & is_neginf(conhist(:, 1:min(nf - 1_IK, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
    end if
    ! Check the values of X, F, CSTRV, CONSTR.
    ! X does not contain NaN if X0 does not and the trust-region/geometry steps are proper.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    ! F cannot be NaN/+Inf due to the moderated extreme barrier.
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    !----------------------------------------------------------------------------------------------!
    ! The following test is not applicable to LINCOA.
    !!if (present(cstrv)) then
    !!    ! CSTRV cannot be NaN/+Inf due to the moderated extreme barrier.
    !!    call assert(.not. (cstrv < 0 .or. is_nan(cstrv) .or. is_posinf(cstrv)), &
    !!        & 'CSTRV is nonnegative and not NaN/+Inf', srname)
    !!end if
    !----------------------------------------------------------------------------------------------!
    if (present(constr)) then
        ! CONSTR cannot contain NaN/-Inf due to the moderated extreme barrier.
        call assert(.not. any(is_nan(constr) .or. is_neginf(constr)), 'CONSTR does not contain NaN/-Inf', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

! Save the history. Note that NF may exceed the maximal amount of history to save. We save X and F
! at the position indexed by MODULO(NF - 1, MAXHIST) + 1. When the solver terminates, the history
! will be reordered so that the information is in the chronological order. Similar for CONSTR, CSTRV.
if (maxxhist > 0) then
    ! We could replace MODULO(NF - 1_IK, MAXXHIST) + 1_IK) with MODULO(NF - 1_IK, MAXHIST) + 1_IK)
    ! based on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
    xhist(:, modulo(nf - 1_IK, maxxhist) + 1_IK) = x
end if
if (maxfhist > 0) then
    fhist(modulo(nf - 1_IK, maxfhist) + 1_IK) = f
end if
if (maxchist > 0) then  ! MAXCHIST > 0 implies PRESENT(CHIST) and PRESENT(CSTRV)
    chist(modulo(nf - 1_IK, maxchist) + 1_IK) = cstrv
end if
if (maxconhist > 0) then  ! MAXCONHIST > 0 implies PRESENT(CONHIST) and PRESENT (CONSTR)
    conhist(:, modulo(nf - 1_IK, maxconhist) + 1_IK) = constr
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then  ! Called after each function evaluation when debugging; can be expensive.
    call assert(size(xhist, 1) == size(x) .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [SIZE(X), MAXXHIST]', srname)
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
        & 'FHIST does not contain NaN/+Inf', srname)
    !----------------------------------------------------------------------------------------------!
    ! The following test is not applicable to LINCOA.
    !!if (present(chist)) then
    !!    call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
    !!    call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) .or. &
    !!        & is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
    !!end if
    !----------------------------------------------------------------------------------------------!
    if (present(conhist) .and. present(constr)) then
        call assert(size(conhist, 1) == size(constr) .and. size(conhist, 2) == maxconhist, &
            & 'SIZE(CONHIST) == [SIZE(CONSTR), MAXCONHIST]', srname)
        call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
            & is_neginf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
    end if
end if

end subroutine savehist


subroutine rangehist(nf, xhist, fhist, chist, conhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine arranges FHIST, XHIST, CHIST, and CONHIST in the chronological order.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
use, non_intrinsic :: debug_mod, only : assert
implicit none

! Inputs
integer(IK), intent(in) :: nf

! In-outputs
real(RP), intent(inout) :: fhist(:)
real(RP), intent(inout) :: xhist(:, :)
real(RP), intent(inout), optional :: chist(:)
real(RP), intent(inout), optional :: conhist(:, :)

! Local variables
integer(IK) :: khist
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: m
integer(IK) :: n
character(len=*), parameter :: srname = 'RANGEHIST'

! Sizes
n = int(size(xhist, 1), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
if (present(chist)) then
    maxchist = int(size(chist), kind(maxchist))
else
    maxchist = 0_IK
end if
if (present(conhist)) then
    m = int(size(conhist, 1), kind(m))
    maxconhist = int(size(conhist, 2), kind(maxconhist))
else
    m = 0_IK
    maxconhist = 0_IK
end if
maxhist = max(maxxhist, maxfhist, maxconhist, maxchist)

! Preconditions
if (DEBUGGING) then
    ! Check the sizes of XHIST, FHIST, CHIST, CONHIST.
    call assert(n >= 1, 'SIZE(XHIST, 1) >= 1', srname)
    call assert(maxxhist * (maxxhist - maxhist) == 0, 'SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
    call assert(maxconhist * (maxconhist - maxhist) == 0, 'SIZE(CONHIST, 2) == 0 or MAXHIST', srname)
    ! Check the values of XHIST, FHIST, CHIST, CONHIST.
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. &
        & is_posinf(fhist(1:min(nf, maxfhist)))), 'FHIST does not contain NaN/+Inf', srname)
    !----------------------------------------------------------------------------------------------!
    ! The following test is not applicable to LINCOA
    !!if (present(chist)) then
    !!call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) .or. &
    !!    & is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
    !!end if
    !----------------------------------------------------------------------------------------------!
    if (present(conhist)) then
        call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
            & is_neginf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

! The ranging should be done only if 0 < MAXXHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxxhist > 0 .and. maxxhist < nf) then
    ! We could replace MODULO(NF - 1_IK, MAXXHIST) + 1_IK) with MODULO(NF - 1_IK, MAXHIST) + 1_IK)
    ! based on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
    khist = modulo(nf - 1_IK, maxxhist) + 1_IK
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
    ! N.B.:
    ! 1. The result of the array constructor is always a rank-1 array (e.g., vector), no matter what
    ! elements are used for the construction.
    ! 2. The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
    ! order of Fortran arrays.
    ! 3. In MATLAB, `xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)]` does the same thing.
end if
! The ranging should be done only if 0 < MAXFHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxfhist > 0 .and. maxfhist < nf) then
    khist = modulo(nf - 1_IK, maxfhist) + 1_IK
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
end if
! The ranging should be done only if 0 < MAXCONHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxconhist > 0 .and. maxconhist < nf) then
    khist = modulo(nf - 1_IK, maxconhist) + 1_IK
    conhist = reshape([conhist(:, khist + 1:maxconhist), conhist(:, 1:khist)], shape(conhist))
end if
! The ranging should be done only if 0 < MAXCHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxchist > 0 .and. maxchist < nf) then
    khist = modulo(nf - 1_IK, maxchist) + 1_IK
    chist = [chist(khist + 1:maxchist), chist(1:khist)]
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
        & 'FHIST does not contain NaN/+Inf', srname)
    if (present(chist)) then
        call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
        !------------------------------------------------------------------------------------------!
        ! The following test is not applicable to LINCOA
        !!call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) .or. &
        !!    & is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
        !------------------------------------------------------------------------------------------!
    end if
    if (present(conhist)) then
        call assert(size(conhist, 1) == m .and. size(conhist, 2) == maxconhist, &
            & 'SIZE(CONHIST) == [M, MAXCONHIST]', srname)
        call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
            & is_neginf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
    end if
end if

end subroutine rangehist


end module history_mod
