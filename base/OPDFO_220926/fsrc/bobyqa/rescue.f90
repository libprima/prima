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
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Monday, September 12, 2022 PM09:22:01
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: rescue


contains


subroutine rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, fhist, fopt, fval,&
    & gopt, hq, pq, sl, su, xbase, xhist, xopt, xpt, zmat, info)
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
! cannot replace any provisional point safely, then set SCORE(KORIG) to -SCORE(KORIG)-SCOREINC;
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
!
! PTSAUX is a 2-by-N real array. For J = 1, 2, ..., N, PTSAUX(1, J) and PTSAUX(2, J) specify the two
! positions of provisional interpolation points when a nonzero step is taken along e_J (the J-th
! coordinate direction) through XBASE + XOPT, as specified below. Usually these steps have length
! DELTA, but other lengths are chosen if necessary in order to satisfy the bound constraints.
!
! PTSID is an integer array of length NPT. Its components denote provisional new positions of the
! interpolation points. The K-th point is a candidate for change if and only if PTSID(K) is nonzero.
! In this case let IP and IQ be the integer parts of PTSID(K) and (PTSID(K)-IP)*(N+1). If IP and IQ
! are both positive, the step from XBASE + XOPT to the new K-th interpolation point is
! PTSAUX(1, IP)*e_IP + PTSAUX(1, IQ)*e_IQ. Otherwise the step is either PTSAUX(1, IP)*e_IP or
! PTSAUX(2, IQ)*e_IQ in the cases IQ=0 or  IP=0, respectively.
!--------------------------------------------------------------------------------------------------!

! Generic modules
use, non_intrinsic :: checkexit_mod, only : checkexit
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, HUGENUM, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: infos_mod, only : MAXFUN_REACHED, INFO_DFT
use, non_intrinsic :: linalg_mod, only : issymmetric, matprod, inprod, r1update, r2update, trueloc
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : hess_mul, setij

! Solver-specif modules
use, non_intrinsic :: update_mod, only : updateh

implicit none

! Inputs
procedure(OBJ) :: calfun
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: maxfun
real(RP), intent(in) :: delta
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: xl(:)  ! XL(N)
real(RP), intent(in) :: xu(:)  ! XU(N)

! In-outputs
integer(IK), intent(inout) :: kopt
integer(IK), intent(inout) :: nf
real(RP), intent(inout) :: bmat(:, :)  !  BMAT(N, NPT + N)
real(RP), intent(inout) :: fhist(:)  ! FHIST(MAXFHIST)
real(RP), intent(inout) :: fopt
real(RP), intent(inout) :: fval(:)  ! FVAL(NPT)
real(RP), intent(inout) :: gopt(:)  ! GOPT(N)
real(RP), intent(inout) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(inout) :: pq(:)  ! PQ(NPT)
real(RP), intent(inout) :: sl(:)  ! SL(N)
real(RP), intent(inout) :: su(:)  ! SU(N)
real(RP), intent(inout) :: xbase(:)  ! XBASE(N)
real(RP), intent(inout) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(inout) :: xopt(:)  ! XOPT(N)
real(RP), intent(inout) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)
! N.B.: BMAT and ZMAT must be INTENT(INOUT) rather than INTENT(OUT); otherwise, they will be
! undefined when the subroutine returns due to NF >= MAXFUN before any calculation starts.

! Outputs
integer(IK), intent(out) :: info

! Local variables
character(len=*), parameter :: srname = 'RESCUE'
integer(IK) :: ip
integer(IK) :: iq
integer(IK) :: ij(2, max(0, size(xpt, 2) - 2 * size(xpt, 1) - 1))
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
real(RP) :: denom
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
real(RP) :: vlag(size(xpt, 1) + size(xpt, 2))
real(RP) :: vlmxsq
real(RP) :: vquad
real(RP) :: wmv(size(xpt, 1) + size(xpt, 2))
real(RP) :: x(size(xpt, 1))
real(RP) :: xp
real(RP) :: v(size(xpt, 1))
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
    call assert(size(xopt) == n, 'SIZE(XOPT) == N, XOPT is finite', srname)
    call assert(all(is_finite(xpt)), 'XPT is finite', srname)
    call assert(all(xpt >= spread(sl, dim=2, ncopies=npt)) .and. &
        & all(xpt <= spread(su, dim=2, ncopies=npt)), 'SL <= XPT <= SU', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NPT+1:NPT+N) is symmetric', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
end if

info = INFO_DFT

! Do nothing if NF already reaches it upper bound.
if (nf >= maxfun) then
    info = MAXFUN_REACHED
    return
end if

! Shift the interpolation points so that XOPT becomes the origin.
xopt = xpt(:, kopt)
xpt = xpt - spread(xopt, dim=2, ncopies=npt)

! Update HQ so that HQ and PQ define the second derivatives of the model after XBASE has been
! shifted to the trust region centre.
v = matprod(xpt, pq) + HALF * sum(pq) * xopt
call r2update(hq, ONE, xopt, v)

! Shift XBASE, SL, SU and XOPT. Set the elements of BMAT and ZMAT to ZERO.
xbase = min(max(xl, xbase + xopt), xu)
sl = min(sl - xopt, ZERO)
su = max(su - xopt, ZERO)
xopt = ZERO

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
! current. Note that, in the BOBYQA paper, the scores are the distances rather than their
! squares. See the paragraph between (5.9) and (5.10) of the BOBYQA paper.
score = sum((xpt)**2, dim=1)
score(kopt) = ZERO  ! Set SCORE(KOPT) to 0 so that KOPT will be skipped when we choose KORIG below.
scoreinc = maxval(score)

! NPROV is the number of provisional points that has not yet been replaced by original points.
nprov = npt - 1_IK

! The following loop runs for at most NPT^2 times: for each value of NPROV, we need at most NPT
! loops to find an original point that can safely replace a provisional point; if such
! a pair of origin and provisional points are found, then NPROV will de reduced by 1; otherwise,
! SCORE will become all zero or negative, and the loop will exit.
do while (any(score(1:npt) > 0) .and. nprov > 0)
    ! Pick the index KORIG of an original point that has not yet replaced one of the provisional
    ! points, giving attention to the closeness to XOPT and to previous tries with KORIG.
    korig = int(minloc(score(1:npt), mask=(score(1:npt) > 0), dim=1), IK)

    ! Calculate VLAG and BETA for the required updating of the H matrix if XPT(:, KORIG) is
    ! reinstated in the set of interpolation points, which means to replace a point in the
    ! following provisional interpolation set XPT_PROV defined in (5.4)--(5.5) of the BOBYQA paper.
    ! 1. XPT_PROV(:, KOPT) = 0;
    ! 2. For each K /= KOPT, if PTSID(K) == 0, then XPT_PROV(:, K) = XPT(:, K); if PTSID(K) > 0,
    ! then XPT_PROV(:, K) has nonzeros only at IP (if IP > 0), IQ (if IQ > 0) positions, where IP
    ! and IQ are the P(J) and Q(J) defined in and below (2.4) of the BOBYQA paper.
    !
    ! First, form the (W - V) vector for XPT(:, KORIG).
    ! In the code below, WMV = W - V = w(XNEW) - w(XOPT) without the (NPT+1)th entry, where
    ! XNEW = XPT(:,KORIG), XOPT = XPT_PROV = 0, and w(.) is defined by (6.3) of the NEWUOA paper
    ! (as well as (4.10) of the BOBYQA paper). Since XOPT = 0, we see from (6.3) that w(XOPT)(K) = 0
    ! for all K except that w(XOPT)(NPT+1) = 1. Thus WMV is the same at w(XNEW) without the
    ! (NPT+1)-th entry. Therefore, WMV= [HALF*MATPROD(XNEW, XPT_PROV)**2, XNEW].
    do k = 1, npt
        if (k == kopt) then
            wmv(k) = ZERO
        else if (ptsid(k) <= 0) then  ! Indeed, PTSID >= 0. So PTSID(K) <= 0 means PTSID(K) = 0.
            wmv(k) = inprod(xpt(:, korig), xpt(:, k))
        else
            ip = int(ptsid(k))  ! IP = 0 if 0 < PTSID(K) < 1. INT(X) = [X], i.e., it rounds X towards 0.
            iq = int(real(n + 1, RP) * ptsid(k) - real(ip * (n + 1), RP))
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
    ! Now calculate VLAG = H*WMV + e_KOPT according to (4.26) of the NEWUOA paper, except VLAG(KOPT).
    vlag(1:npt) = matprod(zmat, matprod(wmv(1:npt), zmat)) + matprod(wmv(npt + 1:npt + n), bmat(:, 1:npt))
    vlag(npt + 1:npt + n) = matprod(bmat, wmv(1:npt + n))
    ! Now calculate BETA. According to (4.12) of the NEWUOA paper (also (4.10) of the BOBYQA paper),
    ! BETA = HALF*|XNEW - XOPT|^4 - WMV'*H*WMV. To calculate WMX'*H*WMV, note that
    ! WMV'*H*WMV = WMV' * [Z*Z', B2^T; B1, B2] * WMV with Z = ZMAT, B1 = BMAT(:, 1:NPT), and
    ! B2 = BMAT(:, NPT+1:NPT+N). Denoting W1 = WMV(1:NPT) and W2 = WMV(NPT+1:NPT+N), we then have
    ! WMV'*H*WMV = |W1*Z|^2 + 2*W1'*B1*W2 + W1'*B2*W2 = |W1*Z|^2 + W1'(B1*W2 + [B1, B2]*WMV).
    bsum = inprod(wmv(1:n), matprod(bmat(:, 1:npt), wmv(1:npt)) + matprod(bmat, wmv))
    beta = HALF * sum(xpt(:, korig)**2)**2 - sum(matprod(wmv(1:npt), zmat)**2) - bsum
    ! Finally, set VLAG(KOPT) to the correct value.
    vlag(kopt) = vlag(kopt) + ONE

    ! For all K with PTSID(K) > 0, calculate the denominator DEN(K) = SIGMA in the updating formula
    ! of H for XPT(:, KORIG) to replace XPT_PROV(:, K).
    hdiag(trueloc(ptsid > 0)) = sum(zmat(trueloc(ptsid > 0), :)**2, dim=2)
    den = ZERO
    den(trueloc(ptsid > 0)) = hdiag(trueloc(ptsid > 0)) * beta + vlag(trueloc(ptsid > 0))**2

    ! KPROV is set to the index of the provisional point that is going to be deleted to make way for
    ! the KORIG-th original point. The choice of KPROV is governed by the avoidance of a small value
    ! of the denominator evaluated above.
    kprov = 0_IK
    denom = ZERO
    if (any(den > 0)) then
        kprov = int(maxloc(den, mask=(.not. is_nan(den)), dim=1), IK)
        denom = den(kprov)
        !!MATLAB: [denom, kprov] = max(den, [], 'omitnan');
    end if
    vlmxsq = HUGENUM
    if (any(.not. is_nan(vlag(1:npt)))) then
        vlmxsq = maxval(vlag(1:npt)**2, mask=(.not. is_nan(vlag(1:npt))))
        !!MATLAB: vlmxsq =  max(vlag(1:npt)**2, [], 'omitnan');
    end if
    if (kopt == 0 .or. denom <= 1.0E-2_RP * vlmxsq) then
        ! Indeed, KOPT == 0 can be removed from the above condition, because KOPT == 0 implies that
        ! DENOM == 0 <= 1.0E-2*VLMXSQ. However, we prefer to mention KOPT == 0 explicitly.
        ! Until finding the next KORIG that renders DENOM > 1.0E-2*VLMXSQ, we will skip the original
        ! interpolation points with a negative or zero score when looking for KORIG (see the
        ! definition of KORIG). When the KORIG satisfying the aforesaid inequality is found, all the
        ! scores will be reset to their absolute values, so that all the original
        ! points with nonzero scores will be considered again.
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
    ! Reset SCORE to ABS(SCORE) so that all the origin points with a nonzero score will be checked
    ! in later loops.
    score = abs(score)

    ! Update the BMAT and ZMAT matrices so that the KORIG-th original point replaces the KORIG-th
    ! provisional point.
    call updateh(korig, beta, vlag, bmat, zmat)

    ! NPROV is the number of provisional points that has not yet been replaced by original points.
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
        if (ptsid(kpt) <= ZERO) then
            cycle
        end if

        ! Absorb PQ(KPT)*XPT(:, KPT)*XPT(:, KPT)^T into the explicit part of the Hessian of the
        ! quadratic model. Implement R1UPDATE properly so that it ensures HQ is symmetric.
        call r1update(hq, pq(kpt), xpt(:, kpt))
        pq(kpt) = ZERO

        ip = int(ptsid(kpt))
        iq = int(real(n + 1, RP) * ptsid(kpt) - real(ip * (n + 1), RP))

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
        x = min(max(xl, xbase + xpt(:, kpt)), xu)
        x(trueloc(xpt(:, kpt) <= sl)) = xl(trueloc(xpt(:, kpt) <= sl))
        x(trueloc(xpt(:, kpt) >= su)) = xu(trueloc(xpt(:, kpt) >= su))
        nf = nf + 1
        call evaluate(calfun, x, f)
        call savehist(nf, x, xhist, f, fhist)
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
            ip = int(ptsid(k))
            iq = int(real(n + 1, RP) * ptsid(k) - real(ip * (n + 1), RP))
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

! Update FOPT, XOPT, and GOPT if necessary.
if (kopt /= kbase) then
    fopt = fval(kopt)
    xopt = xpt(:, kopt)
    gopt = gopt + hess_mul(xopt, xpt, pq, hq)
end if

! Postconditions
if (DEBUGGING) then
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(.not. (is_nan(fopt) .or. is_posinf(fopt)), 'FOPT is not NaN/+Inf', srname)
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
    call assert(size(xopt) == n .and. all(is_finite(xopt)), 'SIZE(XOPT) == N, XOPT is finite', srname)
    call assert(all(xopt >= sl .and. xopt <= su), 'SL <= XOPT <= SU', srname)
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


end module rescue_mod
