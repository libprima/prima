module rescue_mod
!--------------------------------------------------------------------------------------------------!
! This module contains the RESCUE subroutine described in Section 5 of the BOBYQA paper.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
!
! N.B.: According to a test on 20220425, the invocations of RESCUE is rare --- it is never invoked
! on CUTEst unconstrained or bound constrained problems with at most 50 variables unless heavy noise
! is imposed on the function evaluation.
!
! Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Sunday, March 12, 2023 PM07:12:29
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: rescue


contains


subroutine rescue(calfun, solver, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, fhist, fval, &
    & gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat, info)
!--------------------------------------------------------------------------------------------------!
! This subroutine implements "the method of RESCUE" introduced in Section 5 of BOBYQA paper. The
! purpose of this subroutine is to replace a few interpolation points by new points in order to
! improve the geometry of the interpolation set and the conditioning of the interpolation system.
! This is done in the following way.
!
! 1. Define a set of "provisional interpolation points" XPT_PROV around the current XOPT. Similar to
! the construction of the initial interpolation set, XPT_PROV is obtained by perturbing XOPT subject
! to the bound constraints along one or two coordinate directions, the latter taking place only if
! NPT >= 2*N+2. See (5.4)--(5.5) of the BOBYQA paper for details. After defining XPT_PROV, set BMAT
! and ZMAT to represent the H matrix defined in (2.7) of the BOBYQA paper corresponding to XPT_PROV.
! N.B.: In the code, XPT_PROV is not formed explicitly, but represented implicitly by XPT, PTSID,
! and PTSAUX.
! 2. For each "original interpolation point" in XPT, check whether it can replace a point in XPT_PROV
! without damaging the geometry of XPT_PROV, which is indicated by the denominator SIGMA in the
! updating formula of H due to the replacement (see (4.9) of the BOBYQA paper). If yes, update
! XPT_PROV by the performing such an replacement. Continue doing this until all the original points
! get reinstated in XPT_PROV, or we cannot find any original point that can safely replace
! a provisional point. Note the following.
! 2.1. Suppose that the KORIG-th original point is going to replace the KPROV-th provisional point.
! Then we first exchange the KORIG-th and KPROV-th provisional points, and then replace the new
! KORIG-th provisional point by the KORIG-th original point. BMAT and ZMAT are updated accordingly.
! In this way, FVAL(KORIG) is consistent with XPT_PROV(:, KORIG), so that we need not update FVAL.
! 2.1. The KOPT-th original point (i.e., XOPT) is always reinstated in XPT_PROV. This is done by
! replacing XPT_PROV(:, 1) with XOPT. This is the first replacement to perform.
! 2.2. After XOPT is reinstated, the original points are ranked according the scores saved in SCORE.
! SCORE is initialized to the squares of the distances from between the original point and XOPT.
! The KORIG is set to the index of the point with the smallest positive score. If XPT(:, KORIG)
! cannot replace any provisional point safely, then set SCORE(KORIG) to -SCORE(KORIG) - SCOREINC;
! otherwise, set SCORE(KORIG) to 0 and all the other scores to their absolute values. In this way,
! an original point that fails to replace any provisional point during the current attempt will be
! skipped until another original point succeeds in doing so; moreover, the failing point will get
! a lower priority in later attempts. An original point that successfully replaces a provisional
! point will not be tried again due to the zero score.
! 2.3. Once a provisional point is replaced, it will be marked (by setting the corresponding entry
! of PTSID to zero) so that it will not be replaced again.
! 3. When the above procedure finishes, normally most original points are reinstated in XPT_PROV, so
! that XPT_PROV differs from XPT only at very few positions. Set XPT to XPT_PROV, update FVAL at the
! new interpolation points by evaluating F, and then quadratic interpolant [GQ, PQ, HQ] accordingly.
!
! At the end of the subroutine, the elements of BMAT and ZMAT are set in a well-conditioned way to
! the values that are appropriate for the new interpolation points. The elements of GOPT, HQ and PQ
! are also revised to the values that are appropriate to the final quadratic model.
!
! The arguments NF, KOPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL GOPT, HQ, PQ, BMAT, ZMAT, SL
! and, SU have the same meanings as the corresponding arguments of BOBYQB on the entry to RESCUE.
! DELTA is the current trust region radius.
! PTSAUX is a 2-by-N real array. For J = 1, 2, ..., N, PTSAUX(1, J) and PTSAUX(2, J) specify the two
!   positions of provisional interpolation points when a nonzero step is taken along e_J (the J-th
!   coordinate direction) through XBASE + XOPT, as specified below. Usually these steps have length
!   DELTA, but other lengths are chosen if necessary in order to satisfy the bound constraints.
! PTSID is an integer array of length NPT. Its components denote provisional new positions of the
!   interpolation points. The K-th point is a candidate for change if and only if PTSID(K) is
!   nonzero. In this case let IP and IQ be the integer parts of PTSID(K) and (PTSID(K)-IP)*(N+1). If
!   IP and IQ are both positive, the step from XBASE + XOPT to the new K-th interpolation point is
!   PTSAUX(1, IP)*e_IP + PTSAUX(1, IQ)*e_IQ. Otherwise the step is either PTSAUX(1, IP)*e_IP or
!   PTSAUX(2, IQ)*e_IQ in the cases IQ=0 or  IP=0, respectively.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: infos_mod, only : MAXFUN_REACHED, INFO_DFT
use, non_intrinsic :: linalg_mod, only : issymmetric, matprod, inprod, r1update, r2update, trueloc
use, non_intrinsic :: output_mod, only : fmsg
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : hess_mul, setij
use, non_intrinsic :: xinbd_mod, only : xinbd

implicit none

! Inputs
procedure(OBJ) :: calfun
character(len=*), intent(in) :: solver
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: delta
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: xl(:)  ! XL(N)
real(RP), intent(in) :: xu(:)  ! XU(N)

! In-outputs
integer(IK), intent(inout) :: kopt
integer(IK), intent(inout) :: nf
real(RP), intent(inout) :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(inout) :: fval(:)  ! FVAL(NPT)
real(RP), intent(inout) :: gopt(:)  ! GOPT(N)
real(RP), intent(inout) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(inout) :: pq(:)  ! PQ(NPT)
real(RP), intent(inout) :: sl(:)  ! SL(N)
real(RP), intent(inout) :: su(:)  ! SU(N)
real(RP), intent(inout) :: xbase(:)  ! XBASE(N)
real(RP), intent(inout) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(inout) :: xpt(:, :)  ! XPT(N, NPT)

! Outputs
integer(IK), intent(out) :: info
real(RP), intent(out) :: bmat(:, :)  !  BMAT(N, NPT + N)
real(RP), intent(out) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Local variables
character(len=*), parameter :: srname = 'RESCUE'
integer(IK) :: ij(2, max(0, size(xpt, 2) - 2 * size(xpt, 1) - 1))
integer(IK) :: ip
integer(IK) :: iq
integer(IK) :: k
integer(IK) :: kbase
integer(IK) :: korig
integer(IK) :: kprov
integer(IK) :: kpt
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: nprov
integer(IK) :: npt
integer(IK) :: subinfo
logical :: mask(size(xpt, 1))
real(RP) :: beta
real(RP) :: bsum
real(RP) :: den(size(xpt, 2))
real(RP) :: f
real(RP) :: fbase
real(RP) :: hdiag(size(xpt, 2))
real(RP) :: moderr
real(RP) :: pqinc(size(xpt, 2))
real(RP) :: ptsaux(2, size(xpt, 1))
real(RP) :: ptsid(size(xpt, 2))
real(RP) :: score(size(xpt, 2))
real(RP) :: scoreinc
real(RP) :: sfrac
real(RP) :: temp
real(RP) :: v(size(xpt, 1))
real(RP) :: vlag(size(xpt, 1) + size(xpt, 2))
real(RP) :: vquad
real(RP) :: wmv(size(xpt, 1) + size(xpt, 2))
real(RP) :: x(size(xpt, 1))
real(RP) :: xopt(size(xpt, 1))
real(RP) :: xp
real(RP) :: xq
real(RP) :: xxpt(size(xpt, 2))

n = int(size(xpt, 1), kind(n))
npt = int(size(xpt, 2), kind(npt))
maxxhist = int(size(xhist, 2), kind(maxxhist))
maxfhist = int(size(fhist), kind(maxfhist))
maxhist = int(max(maxxhist, maxfhist), kind(maxhist))

! Preconditions
if (DEBUGGING) then
    call assert(abs(iprint) <= 3, 'IPRINT is 0, 1, -1, 2, -2, 3, or -3', srname)
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(maxfun >= npt + 1, 'MAXFUN >= NPT+1', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(delta > 0, 'DELTA > 0', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN/+Inf', srname)
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) is the smallest in FVAL', srname)
    call assert(maxfhist * (maxfhist - maxhist) == 0, 'SIZE(FHIST) == 0 or MAXHIST', srname)
    call assert(size(xl) == n .and. size(xu) == n, 'SIZE(XL) == N == SIZE(XU)', srname)
    call assert(size(sl) == n .and. all(sl <= 0), 'SIZE(SL) == N, SL <= 0', srname)
    call assert(size(su) == n .and. all(su >= 0), 'SIZE(SU) == N, SU >= 0', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) == N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is n-by-n and symmetric', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(all(xbase >= xl .and. xbase <= xu), 'XL <= XBASE <= XU', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(all(xpt >= spread(sl, dim=2, ncopies=npt)) .and. &
        & all(xpt <= spread(su, dim=2, ncopies=npt)), 'SL <= XPT <= SU', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
end if

!====================!
! Calculation starts !
!====================!

info = INFO_DFT

! Do nothing if NF already reaches it upper bound.
! To please Fortran compilers, set BMAT and ZMAT before returning, though they will not be used.
if (nf >= maxfun) then
    bmat = ZERO
    zmat = ZERO
    info = MAXFUN_REACHED
    return
end if

! Shift the interpolation points so that XOPT becomes the origin.
xopt = xpt(:, kopt)
xpt = xpt - spread(xopt, dim=2, ncopies=npt)
xpt(:, kopt) = ZERO

! Update HQ so that HQ and PQ define the second derivatives of the model after XBASE has been
! shifted to the trust region centre.
v = matprod(xpt, pq) + HALF * sum(pq) * xopt
call r2update(hq, ONE, xopt, v)

! Shift XBASE, SL, and SU. Set the elements of BMAT and ZMAT to ZERO.
xbase = min(max(xl, xbase + xopt), xu)
sl = min(xl - xbase, ZERO)
su = max(xu - xbase, ZERO)

! Set the elements of PTSAUX.
ptsaux(1, :) = min(delta, su)
ptsaux(2, :) = max(-delta, sl)
mask = (ptsaux(1, :) + ptsaux(2, :) < 0)
ptsaux([1, 2], trueloc(mask)) = ptsaux([2, 1], trueloc(mask))
mask = (abs(ptsaux(2, :)) < HALF * abs(ptsaux(1, :)))
ptsaux(2, trueloc(mask)) = HALF * ptsaux(1, trueloc(mask))

! Set the identifiers of the artificial interpolation points that are along a coordinate direction
! from XOPT, and set the corresponding nonzero elements of BMAT and ZMAT.
sfrac = HALF / real(n + 1, RP)
ptsid(1) = sfrac
bmat = ZERO
zmat = ZERO
do k = 1, n
    ptsid(k + 1) = real(k, RP) + sfrac
    if (k <= npt - n - 1) then
        ptsid(k + n + 1) = real(k, RP) / real(n + 1, RP) + sfrac
        temp = ONE / (ptsaux(1, k) - ptsaux(2, k))
        bmat(k, k + 1) = -temp + ONE / ptsaux(1, k)
        bmat(k, k + n + 1) = temp + ONE / ptsaux(2, k)
        bmat(k, 1) = -bmat(k, k + 1) - bmat(k, k + n + 1)
        zmat(1, k) = sqrt(TWO) / abs(ptsaux(1, k) * ptsaux(2, k))
        zmat(k + 1, k) = zmat(1, k) * ptsaux(2, k) * temp
        zmat(k + n + 1, k) = -zmat(1, k) * ptsaux(1, k) * temp
    else
        bmat(k, 1) = -ONE / ptsaux(1, k)
        bmat(k, k + 1) = ONE / ptsaux(1, k)
        bmat(k, k + npt) = -HALF * ptsaux(1, k)**2
    end if
end do

! Set any remaining identifiers with their nonzero elements of ZMAT.
ij = setij(n, npt)
do k = 2_IK * n + 2_IK, npt
    ip = ij(1, k - 2 * n - 1)
    iq = ij(2, k - 2 * n - 1)
    ptsid(k) = real(ip, RP) + real(iq, RP) / real(n + 1, RP) + sfrac
    temp = ONE / (ptsaux(1, ip) * ptsaux(1, iq))
    zmat([1_IK, k], k - n - 1) = temp
    zmat([ip + 1, iq + 1], k - n - 1) = -temp
end do

! Update BMAT, ZMAT, ans PTSID so that the 1st and the KOPT-th provisional points are exchanged.
! After the exchanging, the KOPT-th point in the provisional set becomes the zero vector, which is
! exactly the KOPT-th original point (after the shift of XBASE at the beginning of the subroutine).
if (kopt /= 1) then
    bmat(:, [1_IK, kopt]) = bmat(:, [kopt, 1_IK])
    zmat([1_IK, kopt], :) = zmat([kopt, 1_IK], :)
end if
ptsid(1) = ptsid(kopt)
ptsid(kopt) = ZERO

! The squares of the distances from XOPT to the other interpolation points are set at SCORE, which
! will be used to define the index KORIG in the loop below.  Increments of SCOREINC may be added
! later to these scores to balance the consideration of the choice of point that is going to become
! current. Note that, in Powell's BOBYQA code, the initial scores are the squares of the distances,
! but there is no square in the BOBYQA paper (see the paragraph between (5.9) and (5.10) of the
! BOBYQA paper). The latter seem to work better in a test on 20221125.
!score = sum(xpt**2, dim=1)  ! Powell's BOBYQA code
score = sqrt(sum(xpt**2, dim=1))  ! Powell's BOBYQA paper
score(kopt) = ZERO  ! Set SCORE(KOPT) to 0 so that KOPT will be skipped when we choose KORIG below.
scoreinc = maxval(score)

! NPROV is the number of provisional points that has not yet been replaced with original points.
nprov = npt - 1_IK

! The following loop runs for at most NPT^2 times: for each value of NPROV, we need at most NPT
! loops to find an original point that can safely replace a provisional point; if such a pair of
! origin and provisional points are found, then NPROV will de reduced by 1; otherwise, SCORE will
! become all zero or negative, and the loop will exit.
do while (any(score > 0) .and. nprov > 1)   ! Retain at least one provisional point.
! !do while (any(score > 0) .and. nprov > 0)  ! Powell's code. It may not take any provisional point.
! !do while (any(score > 0) .and. nprov > 2)  ! Retain at least two provisional points.
    ! Pick the index KORIG of an original point that has not yet replaced one of the provisional
    ! points, giving attention to the closeness to XOPT and to previous tries with KORIG.
    korig = int(minloc(score, mask=(score > 0), dim=1), kind(korig))

    ! Calculate VLAG and BETA for the required updating of the H matrix if XPT(:, KORIG) is
    ! reinstated in the set of interpolation points, which means to replace a point in the
    ! following provisional interpolation set XPT_PROV defined in (5.4)--(5.5) of the BOBYQA paper.
    ! 1. XPT_PROV(:, KOPT) = 0;
    ! 2. For each K /= KOPT, if PTSID(K) == 0, then XPT_PROV(:, K) = XPT(:, K); if PTSID(K) > 0,
    ! then XPT_PROV(:, K) has nonzeros only at IP (if IP > 0), IQ (if IQ > 0) positions, where IP
    ! and IQ are the P(J) and Q(J) defined in and below (2.4) of the BOBYQA paper.

    ! First, form the (W - V) vector for XPT(:, KORIG).
    ! In the code below, WMV = W - V = w(XNEW) - w(XOPT) without the (NPT+1)th entry, where
    ! XNEW = XPT(:, KORIG), XOPT = XPT_PROV = 0, and w(.) is defined by (6.3) of the NEWUOA paper
    ! (as well as (4.10) of the BOBYQA paper). Since XOPT = 0, we see from (6.3) that w(XOPT)(K) = 0
    ! for all K except that w(XOPT)(NPT+1) = 1. Thus WMV is the same at w(XNEW) without the
    ! (NPT+1)-th entry. Therefore, WMV= [HALF*MATPROD(XNEW, XPT_PROV)**2, XNEW].
    do k = 1, npt
        if (k == kopt) then
            wmv(k) = ZERO
        else if (ptsid(k) <= 0) then  ! Indeed, PTSID >= 0. So PTSID(K) <= 0 means PTSID(K) = 0.
            wmv(k) = inprod(xpt(:, korig), xpt(:, k))
        else
            ip = floor(ptsid(k), kind(ip))  ! IP = 0 if 0 < PTSID(K) < 1.
            iq = floor(real(n + 1, RP) * ptsid(k) - real((n + 1) * ip, RP), kind(iq))
            call assert(ip >= 0 .and. ip <= npt .and. iq >= 0 .and. iq <= npt, '0 <= IP, IQ <= NPT', srname)
            if (ip > 0 .and. iq > 0) then
                wmv(k) = xpt(ip, korig) * ptsaux(1, ip) + xpt(iq, korig) * ptsaux(1, iq)
            elseif (ip > 0) then
                wmv(k) = xpt(ip, korig) * ptsaux(1, ip)
            elseif (iq > 0) then
                wmv(k) = xpt(iq, korig) * ptsaux(2, iq)
            else
                wmv(k) = ZERO
            end if
        end if
        wmv(k) = HALF * wmv(k) * wmv(k)
    end do
    wmv(npt + 1:npt + n) = xpt(:, korig)

    ! Now calculate VLAG = H*WMV + e_KOPT according to (4.26) of the NEWUOA paper except VLAG(KOPT).
    vlag(1:npt) = matprod(zmat, matprod(wmv(1:npt), zmat)) + matprod(wmv(npt + 1:npt + n), bmat(:, 1:npt))
    vlag(npt + 1:npt + n) = matprod(bmat, wmv(1:npt + n))

    ! Now calculate BETA. According to (4.12) of the NEWUOA paper (also (4.10) of the BOBYQA paper),
    ! BETA = HALF*||XNEW - XOPT||^4 - WMV'*H*WMV. To calculate WMX'*H*WMV, note that
    ! WMV'*H*WMV = WMV' * [Z*Z', B2^T; B1, B2] * WMV with Z = ZMAT, B1 = BMAT(:, 1:NPT), and
    ! B2 = BMAT(:, NPT+1:NPT+N). Denoting W1 = WMV(1:NPT) and W2 = WMV(NPT+1:NPT+N), we then have
    ! WMV'*H*WMV = ||W1'*Z||^2 + 2*W1'*B1*W2 + W1'*B2*W2 = ||W1'*Z||^2 + W1'(B1*W2 + [B1, B2]*WMV).
    bsum = inprod(wmv(1:n), matprod(bmat(:, 1:npt), wmv(1:npt)) + matprod(bmat, wmv))
    beta = HALF * sum(xpt(:, korig)**2)**2 - sum(matprod(wmv(1:npt), zmat)**2) - bsum

    ! Finally, set VLAG(KOPT) to the correct value.
    vlag(kopt) = vlag(kopt) + ONE

    ! For all K with PTSID(K) > 0, calculate the denominator DEN(K) = SIGMA in the updating formula
    ! of H for XPT(:, KORIG) to replace XPT_PROV(:, K).
    den = ZERO
    hdiag(trueloc(ptsid > 0)) = sum(zmat(trueloc(ptsid > 0), :)**2, dim=2)
    den(trueloc(ptsid > 0)) = hdiag(trueloc(ptsid > 0)) * beta + vlag(trueloc(ptsid > 0))**2

    ! Attempt setting KPROV to the index of the provisional point to be replaced with the KORIG-th
    ! original interpolation point. We choose KPROV by maximizing DEN(KPROV), which will be the
    ! denominator SIGMA in the updating formula (4.9). In order to avoid a small denominator, we
    ! consider it proper to replace the KPROV-th provisional point with the KORIG-th original point
    ! only if DEN(KPROV) = MAXVAL(DEN) > C*MAXVAL(VLAG(1:NPT)**2), where C is a relatively small
    ! positive constant --- C = 1 is achievable if the rounding errors were not severe; Powell took
    ! C = 0.01, which prefers strongly the original point to the provisional (new) point, as the
    ! latter necessitate new function evaluations. If this inequality is not achievable for the
    ! current KORIG, then we will update SCORE(KORIG) to a negative value and continue the loop with
    ! the next KORIG, which is set to MINLOC(SCORE, MASK=(SCORE > 0)) at the beginning of the loop,
    ! skipping the original interpolation points with a nonpositive score. When a KORIG rendering
    ! the aforesaid inequality is found, SCORE(KORIG) will be set to zero, and all the scores will
    ! be reset to their absolute values, so that future attempts will try the original points that
    ! have not succeeded in replacing a provisional point. The update of SCORE reflects an adaptive
    ! ranking of the original points: points that are closer to XOPT have higher priority, and a
    ! point will be ranked lower if it fails to fulfill MAXVAL(DEN) > C*MAXVAL(VLAG(1:NPT)**2).
    ! Even if KORIG cannot satisfy this condition for now, it may validate the inequality in future
    ! attempts, as BMAT and ZMAT will be updated.
    if (is_finite(sum(abs(vlag))) .and. any(den > 5.0E-2_RP * maxval(vlag(1:npt)**2))) then
        ! The above condition works a bit better than Powell's version below due to the factor 0.05.
        ! !if (any(den > 1.0E-2_RP * maxval(vlag(1:npt)**2))) then  ! Powell' code
        kprov = int(maxloc(den, mask=(.not. is_nan(den)), dim=1), kind(kprov))
        !!MATLAB: [~, kprov] = max(den, [], 'omitnan');
    else
        score(korig) = -score(korig) - scoreinc
        cycle
    end if

    ! Update BMAT, ZMAT, VLAG, and PTSID to exchange the KPROV-th and KORIG-th provisional points.
    ! After the exchanging, the KORIG-th original point will replace the KORIG-th provisional point.
    if (kprov /= korig) then
        bmat(:, [kprov, korig]) = bmat(:, [korig, kprov])
        zmat([kprov, korig], :) = zmat([korig, kprov], :)
        vlag([kprov, korig]) = vlag([korig, kprov])
    end if
    ptsid(kprov) = ptsid(korig)

    ! Set PTSID(KORIG) = 0 so that the KORIG-th provisional point (after the exchanging) will be
    ! skipped in the later loops.
    ptsid(korig) = ZERO
    ! Set SCORE(KORIG) = 0 so that the KORIG-th original point will be skipped in later loops.
    score(korig) = ZERO
    ! Reset SCORE to ABS(SCORE) so that all the original points with a nonzero score will be checked
    ! in later loops.
    score = abs(score)

    ! Update the BMAT and ZMAT matrices so that the KORIG-th original point replaces the KORIG-th
    ! provisional point.
    call updateh_rsc(korig, beta, vlag, bmat, zmat)

    ! NPROV is the number of provisional points that has not yet been replaced with original points.
    nprov = nprov - 1_IK
end do

! All the final positions of the interpolation points have been chosen although any changes have not
! been included yet in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart from the
! shift of XBASE, the updating of the quadratic model remains to be done. The following cycle
! through the new interpolation points begins by putting the new point in XPT(:, KPT) and by setting
! PQ(KPT) to zero. A return occurs if MAXFUN prohibits another value of F or when all the new
! interpolation points are included in the model.
kbase = kopt
fbase = fval(kopt)
if (nprov > 0) then
    do kpt = 1, npt
        if (ptsid(kpt) <= 0) then
            cycle
        end if

        ! Absorb PQ(KPT)*XPT(:, KPT)*XPT(:, KPT)^T into the explicit part of the Hessian of the
        ! quadratic model. Implement R1UPDATE properly so that it ensures HQ is symmetric.
        call r1update(hq, pq(kpt), xpt(:, kpt))
        pq(kpt) = ZERO

        ip = floor(ptsid(kpt), kind(ip))
        iq = floor(real(n + 1, RP) * ptsid(kpt) - real((n + 1) * ip, RP), kind(iq))

        ! Update XPT(:, KPT) to the new point. It contains at most two nonzeros XP and XQ at the IP
        ! and IQ entries.
        xpt(:, kpt) = ZERO
        if (ip > 0 .and. iq > 0) then
            xp = ptsaux(1, ip)
            xpt(ip, kpt) = xp
            xq = ptsaux(1, iq)
            xpt(iq, kpt) = xq
        elseif (ip > 0) then  ! IP > 0, IQ == 0
            xp = ptsaux(1, ip)
            xpt(ip, kpt) = xp
        elseif (iq > 0) then  ! IP == 0, IQ > 0
            xq = ptsaux(2, iq)
            xpt(iq, kpt) = xq
        end if

        ! Calculate F at the new interpolation point, and set MODERR to the factor that is going to
        ! multiply the KPT-th Lagrange function when the model is updated to provide interpolation
        ! to the new function value.
        x = xinbd(xbase, xpt(:, kpt), xl, xu, sl, su)  ! In precise arithmetic, X = XBASE + XPT(:, KPT).
        call evaluate(calfun, x, f)
        nf = nf + 1_IK

        ! Print a message about the function evaluation according to IPRINT.
        call fmsg(solver, iprint, nf, f, x)
        ! Save X, F into the history.
        call savehist(nf, x, xhist, f, fhist)

        ! Update FVAL and KOPT.
        fval(kpt) = f
        if (f < fval(kopt)) then
            kopt = kpt
        end if

        ! Check whether to exit
        subinfo = checkexit(maxfun, nf, f, ftarget, x)
        if (subinfo /= INFO_DFT) then
            info = subinfo
            exit
        end if

        ! Set VQUAD to the value of the current model at the new XPT(:, KPT), which has at most two
        ! nonzeros XP and XQ at the IP and IQ entries respectively.
        vquad = fbase
        if (ip > 0 .and. iq > 0) then
            vquad = vquad + xp * (gopt(ip) + HALF * xp * hq(ip, ip))
            vquad = vquad + xq * (gopt(iq) + HALF * xq * hq(iq, iq))
            vquad = vquad + xp * xq * hq(ip, iq)
            xxpt = xp * xpt(ip, :) + xq * xpt(iq, :)
        elseif (ip > 0) then  ! IP > 0, IQ == 0
            vquad = vquad + xp * (gopt(ip) + HALF * xp * hq(ip, ip))
            xxpt = xp * xpt(ip, :)
        elseif (iq > 0) then  ! IP == 0, IQ > 0
            vquad = vquad + xq * (gopt(iq) + HALF * xq * hq(iq, iq))
            xxpt = xq * xpt(iq, :)
        end if
        vquad = vquad + HALF * inprod(xxpt, pq * xxpt)
        ! N.B.: INPROD(XXPT, PQ * XXPT) = INPROD(X, HESS_MUL(X, XPT, PQ))

        ! Update the quadratic model.
        moderr = f - vquad
        gopt = gopt + moderr * bmat(:, kpt)
        pqinc = moderr * matprod(zmat, zmat(kpt, :))
        pq(trueloc(ptsid <= 0)) = pq(trueloc(ptsid <= 0)) + pqinc(trueloc(ptsid <= 0))
        do k = 1, npt
            if (ptsid(k) <= 0) then
                cycle
            end if
            ip = floor(ptsid(k), kind(ip))
            iq = floor(real(n + 1, RP) * ptsid(k) - real((n + 1) * ip, RP), kind(iq))
            if (ip > 0 .and. iq > 0) then
                hq(ip, ip) = hq(ip, ip) + pqinc(k) * ptsaux(1, ip)**2
                hq(iq, iq) = hq(iq, iq) + pqinc(k) * ptsaux(1, iq)**2
                hq(ip, iq) = hq(ip, iq) + pqinc(k) * ptsaux(1, ip) * ptsaux(1, iq)
                hq(iq, ip) = hq(ip, iq)
            elseif (ip > 0) then  ! IP > 0, IQ == 0
                hq(ip, ip) = hq(ip, ip) + pqinc(k) * ptsaux(1, ip)**2
            elseif (iq > 0) then  ! IP == 0, IP > 0
                hq(iq, iq) = hq(iq, iq) + pqinc(k) * ptsaux(2, iq)**2
            end if
        end do
        ptsid(kpt) = ZERO
    end do
end if

! Update GOPT if necessary.
if (kopt /= kbase) then
    gopt = gopt + hess_mul(xpt(:, kopt), xpt, pq, hq)
end if

!--------------------------------------------------------------------------------------------------!
! Zaikun 20221123: What if we rebuild the model? It seems to worsen the performance of BOBYQA.
! !hq = ZERO
! !pq = omega_mul(1_IK, zmat, fval - fval(kopt))
! !gopt = matprod(bmat(:, 1:npt), fval - fval(kopt)) + hess_mul(xpt(:, kopt), xpt, pq)
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Zaikun 20221123: Shouldn't we correct the models using the new [BMAT, ZMAT]?!
! In this way, we do not even need the quadratic model received by RESCUE is an interpolant.
!real(RP) :: qval(size(xpt, 2))
!qval = [(quadinc(xpt(:, k) - xpt(:, kopt), xpt, gopt, pq, hq), k=1, npt)]
!pq = pq + omega_mul(1_IK, zmat, fval - qval - fval(kopt))
!gopt = gopt + matprod(bmat(:, 1:npt), fval - qval - fval(kopt)) + hess_mul(xpt(:, kopt), xpt, pq)
!--------------------------------------------------------------------------------------------------!

! Postconditions
if (DEBUGGING) then
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(size(fhist) == maxfhist, 'SIZE(FHIST) == MAXFHIST', srname)
    call assert(size(fval) == npt .and. .not. any(is_nan(fval) .or. is_posinf(fval)), &
        & 'SIZE(FVAL) == NPT and FVAL is not NaN/+Inf', srname)
    call assert(.not. any(fval < fval(kopt)), 'FVAL(KOPT) is the smallest in FVAL', srname)
    call assert(size(sl) == n .and. size(su) == n, 'SIZE(SL) == N == SIZE(SU)', srname)
    call assert(size(gopt) == n, 'SIZE(GOPT) == N', srname)
    call assert(size(hq, 1) == n .and. issymmetric(hq), 'HQ is n-by-n and symmetric', srname)
    call assert(size(pq) == npt, 'SIZE(PQ) == NPT', srname)
    call assert(size(xbase) == n .and. all(is_finite(xbase)), 'SIZE(XBASE) == N, XBASE is finite', srname)
    call assert(all(xbase >= xl .and. xbase <= xu), 'XL <= XBASE <= XU', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
    call assert(size(xpt, 1) == n .and. size(xpt, 2) == npt, 'SIZE(XPT) == [N, NPT]', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(all(xpt >= spread(sl, dim=2, ncopies=npt)) .and. &
        & all(xpt <= spread(su, dim=2, ncopies=npt)), 'SL <= XPT <= SU', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(xhist, 1) == n .and. maxxhist * (maxxhist - maxhist) == 0, &
        & 'SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST', srname)
end if

end subroutine rescue


subroutine updateh_rsc(knew, beta, vlag_in, bmat, zmat, info)
! !!! N.B.: UPDATEH_RSC is only used by RESCUE.
!--------------------------------------------------------------------------------------------------!
! This subroutine updates arrays BMAT and ZMAT in order to replace the interpolation point
! XPT(:, KNEW) by XNEW = XPT(:, KOPT) + D. See Section 4 of the BOBYQA paper. [BMAT, ZMAT] describes
! the matrix H in the BOBYQA paper (eq. 2.7), which is the inverse of the coefficient matrix of the
! KKT system for the least-Frobenius norm interpolation problem: ZMAT holds a factorization of the
! leading NPT*NPT submatrix OMEGA of H, the factorization being OMEGA = ZMAT*ZMAT^T; BMAT holds the
! last N ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row and (NPT + 1)th
! column of H are not stored as they are unnecessary for the calculation.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: infnan_mod, only : is_finite
use, non_intrinsic :: infos_mod, only : INFO_DFT, DAMAGING_ROUNDING
use, non_intrinsic :: linalg_mod, only : planerot, matprod, outprod, symmetrize, issymmetric
implicit none

! Inputs
integer(IK), intent(in) :: knew
real(RP), intent(in) :: beta
real(RP), intent(in) :: vlag_in(:)  ! VLAG(NPT + N)

! In-outputs
real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Outputs
integer(IK), intent(out), optional :: info

! Local variables
character(len=*), parameter :: srname = 'UPDATEH_RSC'
integer(IK) :: j
integer(IK) :: n
integer(IK) :: npt
real(RP) :: alpha
real(RP) :: denom
real(RP) :: grot(2, 2)
real(RP) :: hcol(size(bmat, 2))
real(RP) :: sqrtdn
real(RP) :: tau
real(RP) :: v1(size(bmat, 1))
real(RP) :: v2(size(bmat, 1))
real(RP) :: vlag(size(vlag_in))
real(RP) :: ztest

! Sizes.
n = int(size(bmat, 1), kind(n))
npt = int(size(bmat, 2) - size(bmat, 1), kind(npt))

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(knew >= 1 .and. knew <= npt, '1 <= KNEW <= NPT', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(vlag_in) == npt + n, 'SIZE(VLAG) == NPT + N', srname)

    ! The following is too expensive to check.
    !tol = 1.0E-2_RP
    !call wassert(errh(bmat, zmat, xpt) <= tol .or. RP == kind(0.0), &
    !    & 'H = W^{-1} in (2.7) of the BOBYQA paper', srname)
end if

!====================!
! Calculation starts !
!====================!

if (present(info)) then
    info = INFO_DFT
end if

! We must not do anything if KNEW is 0. This can only happen sometimes after a trust-region step.
if (knew <= 0) then  ! KNEW < 0 is impossible if the input is correct.
    return
end if

! Read VLAG, and calculate parameters for the updating formula (4.9) and (4.14) of the BOBYQA paper.
vlag = vlag_in
tau = vlag(knew)
! In theory, DENOM can also be calculated after ZMAT is rotated below. However, this worsened the
! performance of BOBYQA in a test on 20220413.
denom = sum(zmat(knew, :)**2) * beta + tau**2

! Quite rarely, due to rounding errors, VLAG or BETA may not be finite, or DENOM may not be
! positive. In such cases, [BMAT, ZMAT] would be destroyed by the update, and hence we would rather
! not update them at all. Or should we simply terminate the algorithm?
if (.not. (is_finite(sum(abs(vlag)) + abs(beta)) .and. denom > 0)) then
    if (present(info)) then
        info = DAMAGING_ROUNDING
    end if
    return
end if

! After the following line, VLAG = H*w - e_KNEW in the NEWUOA paper (where t = KNEW).
vlag(knew) = vlag(knew) - ONE

! Apply Givens rotations to put zeros in the KNEW-th row of ZMAT. After this, ZMAT(KNEW, :) contains
! only one nonzero at ZMAT(KNEW, 1). Entries of ZMAT are treated as 0 if the moduli are at most ZTEST.
ztest = 1.0E-20_RP * maxval(abs(zmat))
do j = 2, npt - n - 1_IK
    if (abs(zmat(knew, j)) > ztest) then
        grot = planerot(zmat(knew, [1_IK, j]))
        zmat(:, [1_IK, j]) = matprod(zmat(:, [1_IK, j]), transpose(grot))
    end if
    zmat(knew, j) = ZERO
end do

! Put the KNEW-th column of the unupdated H (except for the (NPT+1)th entry) into HCOL.
hcol(1:npt) = zmat(knew, 1) * zmat(:, 1)
hcol(npt + 1:npt + n) = bmat(:, knew)

! Complete the updating of ZMAT. See (4.14) of the BOBYQA paper.
sqrtdn = sqrt(denom)
zmat(:, 1) = (tau / sqrtdn) * zmat(:, 1) - (zmat(knew, 1) / sqrtdn) * vlag(1:npt)

! Finally, update the matrix BMAT. It implements the last N rows of (4.9) in the BOBYQA paper.
alpha = hcol(knew)
v1 = (alpha * vlag(npt + 1:npt + n) - tau * hcol(npt + 1:npt + n)) / denom
v2 = (-beta * hcol(npt + 1:npt + n) - tau * vlag(npt + 1:npt + n)) / denom
bmat = bmat + outprod(v1, vlag) + outprod(v2, hcol) !call r2update(bmat, ONE, v1, vlag, ONE, v2, hcol)
! Numerically, the update above does not guarantee BMAT(:, NPT+1 : NPT+N) to be symmetric.
call symmetrize(bmat(:, npt + 1:npt + n))

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT)==[N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)

    ! The following is too expensive to check.
    ! !if (n * npt <= 50) then
    ! !    xpt_test = xpt
    ! !    xpt_test(:, knew) = xpt(:, kopt) + d
    ! !    call assert(errh(bmat, zmat, xpt_test) <= tol .or. RP == kind(0.0), &
    ! !        & 'H = W^{-1} in (2.7) of the BOBYQA paper', srname)
    ! !end if
end if
end subroutine updateh_rsc


end module rescue_mod
