module history_mod
!--------------------------------------------------------------------------------------------------!
! This module provides subroutines that handle the X/F/C histories of the solver, taking into
! account that MAXHIST may be smaller than NF.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Friday, December 17, 2021 PM11:55:58
!--------------------------------------------------------------------------------------------------!

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
