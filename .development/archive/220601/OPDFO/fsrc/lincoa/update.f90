module update_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the update of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Thursday, May 26, 2022 AM11:10:19
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: update


contains


subroutine update(kopt, step, xpt, idz, knew, bmat, zmat)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, wassert
use, non_intrinsic :: infnan_mod, only : is_nan, is_finite
use, non_intrinsic :: powalg_mod, only : updateh, calden

implicit none

! Inputs
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: step(:)  ! STEP(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! In-outputs
integer(IK), intent(inout) :: idz
integer(IK), intent(inout) :: knew
real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Local variables
character(len=*), parameter :: srname = 'UPDATE'
real(RP) :: denabs(size(xpt, 2)), distsq(size(xpt, 2)), score(size(xpt, 2))
real(RP) :: weight(size(xpt, 2))
integer(IK) :: n
integer(IK) :: npt
!real(RP) :: xopt(size(xpt, 1))
!real(RP) :: xdist(size(xpt, 2))!, xxpt(size(xpt, 2))!, sxpt(size(xpt, 2)), vtmp(size(xpt, 2))


! Sizes.
n = size(xpt, 1)
npt = size(xpt, 2)

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call wassert(size(step) == n .and. all(is_finite(step)), 'SIZE(STEP) == N, STEP is finite', srname)
    call wassert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
end if

!
!     The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM, XSXPT and STEP are
!       identical to the corresponding arguments in SUBROUTINE LINCOB.
!     KOPT is such that XPT(KOPT,.) is the current trust region centre.
!     KNEW on exit is usually positive, and then it is the index of an
!       interpolation point to be moved to the position XPT(KOPT,.)+STEP(.).
!       It is set on entry either to its final value or to 0. In the latter
!       case, the final value of KNEW is chosen to maximize the denominator
!       of the matrix updating formula times a weighting factor.
!     VLAG and W are used for working space, the first NPT+N elements of
!       both of these vectors being required.
!
!     The arrays BMAT and ZMAT with IDZ are updated, the new matrices being
!       the ones that are suitable after the shift of the KNEW-th point to
!       the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to ZERO
!       occurs if the calculation fails due to a ZERO denominator in the
!       updating formula, which should never happen.
!
!     Set some constants.
!
!
!     Calculate VLAG and BETA for the current choice of STEP. The first NPT
!       elements of VLAG are set to the values of the Lagrange functions at
!       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
!       in W, where W_check is defined in a paper on the updating method.
!


! If KNEW is ZERO initially, then pick the index of the interpolation point to be deleted, by
! maximizing the absolute value of the denominator of the updating formula times a weighting factor.
if (knew == 0) then
    denabs = abs(calden(kopt, bmat, step, xpt, zmat, idz))
    distsq = sum((xpt - spread(xpt(:, kopt), dim=2, ncopies=npt))**2, dim=1)
    weight = distsq**2
    score = weight * denabs

    if (any(score > 0)) then
        ! SCORE(K) is NaN implies DENABS(K) is NaN, but we want DENABS to be big. So we exclude such K.
        knew = int(maxloc(score, mask=(.not. is_nan(score)), dim=1), IK)
        !!MATLAB: [~, knew] = max(score, [], 'omitnan');
    else
        ! Powell's code does not handle this case, leaving KNEW = 0 and leading to a segfault.
        knew = int(maxloc(distsq, dim=1), IK)
    end if
end if
call assert(1 <= knew .and. knew <= npt, '1 <= KNEW <= NPT', srname)

call updateh(knew, kopt, idz, step, xpt, bmat, zmat)
end subroutine update


end module update_mod
