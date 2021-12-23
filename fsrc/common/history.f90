module history_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines that handle the X/F/C histories of the solver, taking into
! account that MAXHIST may be smaller than NF.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Thursday, December 23, 2021 AM11:05:09
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: prehist, savehist, rangehist

interface prehist
    module procedure prehist_unc, prehist_nlc
end interface prehist

interface savehist
    module procedure savehist_unc, savehist_nlc
end interface savehist

interface rangehist
    module procedure rangehist_unc, rangehist_nlc
end interface rangehist


contains


subroutine prehist_unc(maxhist, n, output_fhist, fhist, output_xhist, xhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine revises MAXHIST according to MAXMEMORY, and allocates memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN) if they are requested; replace MAXFUN with 0 for
! the history that is not requested.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, MAXMEMORY, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : int
use, non_intrinsic :: memory_mod, only : safealloc, cstyle_sizeof
implicit none

! Inputs
integer(IK), intent(in) :: n
logical, intent(in) :: output_fhist
logical, intent(in) :: output_xhist

! In-outputs
integer(IK), intent(inout) :: maxhist

! Outputs
real(RP), intent(out), allocatable :: fhist(:)
real(RP), intent(out), allocatable :: xhist(:, :)

! Local variables
character(len=*), parameter :: srname = 'PREHIST_UNC'
integer(IK) :: maxhist_in
integer(IK) :: unit_memo

! Preconditions
if (DEBUGGING) then
    call assert(maxhist >= 0, 'MAXHIST >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
end if

!====================!
! Calculation starts !
!====================!

! Save the input value of MAXHIST for debugging.
maxhist_in = maxhist

! Revise MAXHIST according to MAXMEMORY, i.e., the maximal memory allowed for the history.
unit_memo = int(output_xhist) * n + int(output_fhist)
unit_memo = unit_memo * cstyle_sizeof(0.0_RP)
if (unit_memo <= 0) then  ! No output of history is requested
    maxhist = 0_IK
elseif (maxhist > MAXMEMORY / unit_memo) then
    maxhist = int(MAXMEMORY / unit_memo, kind(maxhist))
    ! We cannot simply set MAXHIST = MIN(MAXHIST, MAXMEMORY/UNIT_MEMO), as they may not have
    ! the same kind, and compilers may complain. We may convert them, but overflow may occur.
end if

call safealloc(fhist, maxhist * int(output_fhist))
call safealloc(xhist, n, maxhist * int(output_xhist))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(maxhist >= 0 .and. maxhist <= maxhist_in, '0 <= MAXHIST <= MAXHIST_IN', srname)
    call assert(maxhist == max(size(fhist), size(xhist, 2)), &
        & 'MAXHIST == MAX(SIZE(CHIST), SIZE(CONHIST, 2), SIZE(FHIST), SIZE(XHIST, 2))', srname)
    call assert(int(maxhist, kind(MAXMEMORY)) * int(unit_memo, kind(MAXMEMORY)) <= MAXMEMORY, &
        & 'the history will not take more memory than MAXMEMORY', srname)
    call assert(allocated(fhist) .and. allocated(xhist), 'the history is allocated', srname)
    call assert(size(fhist) == maxhist * int(output_fhist), &
        & 'if FHIST is requested, then SIZE(FHIST) == MAXHIST; otherwise, SIZE(FHIST) == 0', srname)
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxhist * int(output_xhist), &
        & 'if XHIST is requested, then SIZE(XHIST) == [N, MAXHIST]; otherwise, SIZE(XHIST) == [N, 0]', srname)
    call assert(int(size(fhist) + size(xhist), kind(MAXMEMORY)) * &
       & int(cstyle_sizeof(0.0_RP), kind(MAXMEMORY)) <= MAXMEMORY, &
        & 'the history will not take more memory than MAXMEMORY', srname)
end if
end subroutine prehist_unc

subroutine prehist_nlc(maxhist, m, n, output_chist, chist, output_conhist, conhist, output_fhist, fhist, output_xhist, xhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine revises MAXHIST according to MAXMEMORY, and allocates memory for the history.
! In MATLAB/Python/Julia/R implementation, we should simply set MAXHIST = MAXFUN and initialize
! CHIST = NaN(1, MAXFUN), CONHIST = NaN(M, MAXFUN), FHIST = NaN(1, MAXFUN), XHIST = NaN(N, MAXFUN)
! if they are requested; replace MAXFUN with 0 for the history that is not requested.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, MAXMEMORY, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : int
use, non_intrinsic :: memory_mod, only : safealloc, cstyle_sizeof
implicit none

! Inputs
integer(IK), intent(in) :: m
integer(IK), intent(in) :: n
logical, intent(in) :: output_chist
logical, intent(in) :: output_conhist
logical, intent(in) :: output_fhist
logical, intent(in) :: output_xhist

! In-outputs
integer(IK), intent(inout) :: maxhist

! Outputs
real(RP), intent(out), allocatable :: chist(:)
real(RP), intent(out), allocatable :: conhist(:, :)
real(RP), intent(out), allocatable :: fhist(:)
real(RP), intent(out), allocatable :: xhist(:, :)

! Local variables
character(len=*), parameter :: srname = 'PREHIST_NLC'
integer(IK) :: maxhist_in
integer(IK) :: unit_memo

! Preconditions
if (DEBUGGING) then
    call assert(maxhist >= 0, 'MAXHIST >= 0', srname)
    call assert(m >= 0, 'M >= 0', srname)
    call assert(n >= 1, 'N >= 1', srname)
end if

!====================!
! Calculation starts !
!====================!

! Save the input value of MAXHIST for debugging.
maxhist_in = maxhist

! Revise MAXHIST according to MAXMEMORY, i.e., the maximal memory allowed for the history.
unit_memo = int(output_xhist) * n + int(output_fhist) + int(output_conhist) * m + int(output_chist)
unit_memo = unit_memo * cstyle_sizeof(0.0_RP)
if (unit_memo <= 0) then  ! No output of history is requested
    maxhist = 0_IK
elseif (maxhist > MAXMEMORY / unit_memo) then
    maxhist = int(MAXMEMORY / unit_memo, kind(maxhist))
    ! We cannot simply set MAXHIST = MIN(MAXHIST, MAXMEMORY/UNIT_MEMO), as they may not have
    ! the same kind, and compilers may complain. We may convert them, but overflow may occur.
end if

call safealloc(chist, maxhist * int(output_chist))
call safealloc(conhist, m, maxhist * int(output_conhist))
call safealloc(fhist, maxhist * int(output_fhist))
call safealloc(xhist, n, maxhist * int(output_xhist))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(maxhist >= 0 .and. maxhist <= maxhist_in, '0 <= MAXHIST <= MAXHIST_IN', srname)
    call assert(maxhist == max(size(chist), size(conhist, 2), size(fhist), size(xhist, 2)), &
        & 'MAXHIST == MAX(SIZE(CHIST), SIZE(CONHIST, 2), SIZE(FHIST), SIZE(XHIST, 2))', srname)
    call assert(int(maxhist, kind(MAXMEMORY)) * int(unit_memo, kind(MAXMEMORY)) <= MAXMEMORY, &
        & 'the history will not take more memory than MAXMEMORY', srname)
    call assert(allocated(chist) .and. allocated(conhist) .and. allocated(fhist) .and. allocated(xhist), &
        & 'the history is allocated', srname)
    call assert(size(chist) == maxhist * int(output_chist), &
        & 'if CHIST is requested, then SIZE(CHIST) == MAXHIST; otherwise, SIZE(CHIST) == 0', srname)
    call assert(size(conhist, 1) == m .and. size(conhist, 2) == maxhist * int(output_conhist), &
        & 'if CONHIST is requested, then SIZE(CONHIST) == [M, MAXHIST]; otherwise, SIZE(CONHIST) == [M, 0]', srname)
    call assert(size(fhist) == maxhist * int(output_fhist), &
        & 'if FHIST is requested, then SIZE(FHIST) == MAXHIST; otherwise, SIZE(FHIST) == 0', srname)
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxhist * int(output_xhist), &
        & 'if XHIST is requested, then SIZE(XHIST) == [N, MAXHIST]; otherwise, SIZE(XHIST) == [N, 0]', srname)
    call assert(int(size(chist) + size(conhist) + size(fhist) + size(xhist), kind(MAXMEMORY)) * &
       & int(cstyle_sizeof(0.0_RP), kind(MAXMEMORY)) <= MAXMEMORY, &
        & 'the history will not take more memory than MAXMEMORY', srname)
end if
end subroutine prehist_nlc


subroutine savehist_unc(nf, f, x, fhist, xhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine saves X and F into XHIST and FHIST respectively.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf
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

! Sizes
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    ! Check the size of X.
    call assert(n >= 1, 'N >= 1', srname)
    ! Check the sizes of XHIST, FHIST.
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    ! Check the values of XHIST, FHIST, up to the (NF - 1)th position.
    ! As long as this subroutine is called, XHIST contain only finite values.
    call assert(all(is_finite(xhist(:, 1:min(nf - 1_IK, maxxhist)))), 'XHIST is finite', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf - 1_IK, maxfhist))) .or. &
        & is_posinf(fhist(1:min(nf - 1_IK, maxfhist)))), 'FHIST does not contain NaN/+Inf', srname)
    ! Check the values of X, F.
    ! X does not contain NaN if X0 does not and the trust-region/geometry steps are proper.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    ! F cannot be NaN/+Inf due to the moderated extreme barrier.
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

! Save the history. Note that NF may exceed the maximal amount of history to save. We save X and F
! at the position indexed by MOD(NF - 1, MAXHIST) + 1. When the solver terminates, the history
! will be re-ordered so that the information is in the chronological order.
if (maxxhist > 0) then
    ! We could replace MOD(NF - 1_IK, MAXXHIST) + 1_IK) with MOD(NF - 1_IK, MAXHIST) + 1_IK) based
    ! on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
    xhist(:, mod(nf - 1_IK, maxxhist) + 1_IK) = x
end if
if (maxfhist > 0) then
    fhist(mod(nf - 1_IK, maxfhist) + 1_IK) = f
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
end if

end subroutine savehist_unc


subroutine savehist_nlc(nf, constr, cstrv, f, x, chist, conhist, fhist, xhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine saves X, F, CONSTR, and CSTRV into XHIST, FHIST, CONHIST, and CHIST respectively.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite, is_posinf, is_neginf
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

! Sizes.
m = int(size(constr), kind(m))
n = int(size(x), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxconhist = int(size(conhist, 2), kind(maxconhist))
maxchist = int(size(chist), kind(maxchist))
maxhist = max(maxxhist, maxfhist, maxconhist, maxchist)

! Preconditions
if (DEBUGGING) then  ! Called after each function evaluation when debugging; can be expensive.
    ! Check the size of X.
    call assert(n >= 1, 'N >= 1', srname)
    ! Check the sizes of XHIST, FHIST, CONHIST, CHIST.
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(conhist, 1) == m .and. maxconhist * (maxconhist - maxhist) == 0, &
        & 'SIZE(CONHIST, 1) == N, SIZE(CONHIST, 2) == 0 or MAXNHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
    ! Check the values of XHIST, FHIST, CONHIST, CHIST, up to the (NF - 1)th position.
    ! As long as this subroutine is called, XHIST contain only finite values.
    call assert(all(is_finite(xhist(:, 1:min(nf - 1_IK, maxxhist)))), 'XHIST is finite', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf - 1_IK, maxfhist))) .or. &
        & is_posinf(fhist(1:min(nf - 1_IK, maxfhist)))), 'FHIST does not contain NaN/+Inf - 1_IK', srname)
    call assert(.not. any(is_nan(conhist(:, 1:min(nf - 1_IK, maxconhist))) .or. &
        & is_neginf(conhist(:, 1:min(nf - 1_IK, maxconhist)))), 'CONHIST does not contain NaN/-Inf - 1_IK', srname)
    call assert(.not. any(chist(1:min(nf - 1_IK, maxchist)) < 0 .or. &
        & is_nan(chist(1:min(nf - 1_IK, maxchist))) .or. is_posinf(chist(1:min(nf - 1_IK, maxchist)))), &
        & 'CHIST does not contain nonnegative values or NaN/+Inf - 1_IK', srname)
    ! Check the values of X, F, CONSTR, CSTRV.
    ! X does not contain NaN if X0 does not and the trust-region/geometry steps are proper.
    call assert(.not. any(is_nan(x)), 'X does not contain NaN', srname)
    ! F cannot be NaN/+Inf due to the moderated extreme barrier.
    call assert(.not. (is_nan(f) .or. is_posinf(f)), 'F is not NaN/+Inf', srname)
    ! CONSTR cannot contain NaN/-Inf due to the moderated extreme barrier.
    call assert(.not. any(is_nan(constr) .or. is_neginf(constr)), 'CONSTR does not contain NaN/-Inf', srname)
    ! CSTRV cannot be NaN/+Inf due to the moderated extreme barrier.
    call assert(.not. (cstrv < 0 .or. is_nan(cstrv) .or. is_posinf(cstrv)), 'CSTRV is nonnegative and not NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

! Save the history. Note that NF may exceed the maximal amount of history to save. We save X, F,
! CONSTR, and CSTRV at the position indexed by MOD(NF - 1, MAXHIST) + 1. When the solver terminates,
! the history will be re-ordered so that the information is in the chronological order.
if (maxxhist > 0) then
    ! We could replace MOD(NF - 1_IK, MAXXHIST) + 1_IK) with MOD(NF - 1_IK, MAXHIST) + 1_IK) based
    ! on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
    xhist(:, mod(nf - 1_IK, maxxhist) + 1_IK) = x
end if
if (maxfhist > 0) then
    fhist(mod(nf - 1_IK, maxfhist) + 1_IK) = f
end if
if (maxconhist > 0) then
    conhist(:, mod(nf - 1_IK, maxconhist) + 1_IK) = constr
end if
if (maxchist > 0) then
    chist(mod(nf - 1_IK, maxchist) + 1_IK) = cstrv
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then  ! Called after each function evaluation when debugging; can be expensive.
    call assert(size(xhist, 1) == n .and. size(xhist, 2) == maxxhist, 'SIZE(XHIST) == [N, MAXXHIST]', srname)
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. is_posinf(fhist(1:min(nf, maxfhist)))), &
        & 'FHIST does not contain NaN/+Inf', srname)
    call assert(size(conhist, 1) == m .and. size(conhist, 2) == maxconhist, &
        & 'SIZE(CONHIST) == [M, MAXCONHIST]', srname)
    call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
        & is_neginf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
    call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
    call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) .or. &
        & is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
end if

end subroutine savehist_nlc


subroutine rangehist_unc(nf, fhist, xhist)
!--------------------------------------------------------------------------------------------------!
! This subroutine arranges FHIST and XHIST in the chronological order.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf
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
integer(IK) :: n
character(len=*), parameter :: srname = 'RANGEHIST_UNC'

! Sizes
n = int(size(xhist, 1), kind(n))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    ! Check the sizes of XHIST, FHIST.
    call assert(n >= 1, 'N >= 1', srname)
    call assert(maxxhist * (maxxhist - maxhist) == 0, 'SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    ! Check the values of XHIST, FHIST.
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. &
        & is_posinf(fhist(1:min(nf, maxfhist)))), 'FHIST does not contain NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

! The ranging should be done only if 0 < MAXXHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxxhist > 0 .and. maxxhist < nf) then
    ! We could replace MOD(NF - 1_IK, MAXXHIST) + 1_IK) with MOD(NF - 1_IK, MAXHIST) + 1_IK) based
    ! on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
    khist = mod(nf - 1_IK, maxxhist) + 1_IK
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
    ! N.B.:
    ! 1. The result of the array constructor is always a rank-1 array (e.g., vector), no matter what
    ! elements are used to construct the array.
    ! 2. The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
    ! order of Fortran arrays.
    ! 3. In MATLAB, `xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)]` does the same thing.
end if
! The ranging should be done only if 0 < MAXFHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxfhist > 0 .and. maxfhist < nf) then
    khist = mod(nf - 1_IK, maxfhist) + 1_IK
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
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
end if

end subroutine rangehist_unc


subroutine rangehist_nlc(nf, chist, conhist, fhist, xhist)
! This subroutine arranges FHIST and XHIST, CONHIST, and CHIST in the chronological order.
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_neginf
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
integer(IK) :: m
integer(IK) :: maxchist
integer(IK) :: maxconhist
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
character(len=*), parameter :: srname = 'RANGEHIST_NLC'

! Sizes
n = int(size(xhist, 1), kind(n))
m = int(size(conhist, 1), kind(m))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxconhist = int(size(conhist, 2), kind(maxconhist))
maxchist = int(size(chist), kind(maxchist))
maxhist = max(maxxhist, maxfhist, maxconhist, maxchist)

! Preconditions
if (DEBUGGING) then
    ! Check the sizes of XHIST, FHIST, CONHIST, CHIST.
    call assert(n >= 1, 'N >= 1', srname)
    call assert(maxxhist * (maxxhist - maxhist) == 0, 'SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(maxconhist * (maxconhist - maxhist) == 0, 'SIZE(CONHIST, 2) == 0 or MAXHIST', srname)
    call assert(maxchist * (maxchist - maxhist) == 0, 'SIZE(CHIST) == 0 or MAXHIST', srname)
    ! Check the values of XHIST, FHIST, CONHIST, CHIST.
    call assert(.not. any(is_nan(xhist(:, 1:min(nf, maxxhist)))), 'XHIST does not contain NaN', srname)
    ! The last calculated X can be Inf (finite + finite can be Inf numerically).
    call assert(.not. any(is_nan(fhist(1:min(nf, maxfhist))) .or. &
        & is_posinf(fhist(1:min(nf, maxfhist)))), 'FHIST does not contain NaN/+Inf', srname)
    call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
        & is_neginf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
    call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) .or. &
        & is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
end if

!====================!
! Calculation starts !
!====================!

! The ranging should be done only if 0 < MAXXHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxxhist > 0 .and. maxxhist < nf) then
    ! We could replace MOD(NF - 1_IK, MAXXHIST) + 1_IK) with MOD(NF - 1_IK, MAXHIST) + 1_IK) based
    ! on the assumption that MAXXHIST == 0 or MAXHIST. For robustness, we do not do that.
    khist = mod(nf - 1_IK, maxxhist) + 1_IK
    xhist = reshape([xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)], shape(xhist))
    ! N.B.:
    ! 1. The result of the array constructor is always a rank-1 array (e.g., vector), no matter what
    ! elements are used to construct the array.
    ! 2. The above combination of SHAPE and RESHAPE fulfills our desire thanks to the COLUMN-MAJOR
    ! order of Fortran arrays.
    ! 3. In MATLAB, `xhist = [xhist(:, khist + 1:maxxhist), xhist(:, 1:khist)]` does the same thing.
end if
! The ranging should be done only if 0 < MAXFHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxfhist > 0 .and. maxfhist < nf) then
    khist = mod(nf - 1_IK, maxfhist) + 1_IK
    fhist = [fhist(khist + 1:maxfhist), fhist(1:khist)]
end if
! The ranging should be done only if 0 < MAXCONHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxconhist > 0 .and. maxconhist < nf) then
    khist = mod(nf - 1_IK, maxconhist) + 1_IK
    conhist = reshape([conhist(:, khist + 1:maxconhist), conhist(:, 1:khist)], shape(conhist))
end if
! The ranging should be done only if 0 < MAXCHIST < NF. Otherwise, it leads to errors/wrong results.
if (maxchist > 0 .and. maxchist < nf) then
    khist = mod(nf - 1_IK, maxchist) + 1_IK
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
    call assert(size(conhist, 1) == m .and. size(conhist, 2) == maxconhist, &
        & 'SIZE(CONHIST) == [M, MAXCONHIST]', srname)
    call assert(.not. any(is_nan(conhist(:, 1:min(nf, maxconhist))) .or. &
        & is_neginf(conhist(:, 1:min(nf, maxconhist)))), 'CONHIST does not contain NaN/-Inf', srname)
    call assert(size(chist) == maxchist, 'SIZE(CHIST) == MAXCHIST', srname)
    call assert(.not. any(chist(1:min(nf, maxchist)) < 0 .or. is_nan(chist(1:min(nf, maxchist))) .or. &
        & is_posinf(chist(1:min(nf, maxchist)))), 'CHIST does not contain nonnegative values or NaN/+Inf', srname)
end if

end subroutine rangehist_nlc


end module history_mod
