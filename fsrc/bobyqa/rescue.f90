module rescue_mod
!--------------------------------------------------------------------------------------------------!
! This module contains the `rescue` subroutine.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the BOBYQA paper.
!
! N.B.: According to a test on 20220425, the invocations of RESCUE is rare --- it is never invoked
! on CUTEst unconstrained or bound constrained problems with at most 50 variables unless heavy noise
! is imposed on the function evaluation.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Wednesday, May 25, 2022 PM10:15:37
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: rescue


contains


subroutine rescue(calfun, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, bmat, fhist, fval, &
    & gopt, hq, pq, sl, su, vlag, xbase, xhist, xopt, xpt, zmat, f)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: evaluate_mod, only : evaluate
use, non_intrinsic :: history_mod, only : savehist
use, non_intrinsic :: infnan_mod, only : is_nan, is_posinf, is_finite
use, non_intrinsic :: linalg_mod, only : issymmetric, matprod, inprod, r1update, r2update, trueloc!, norm
use, non_intrinsic :: pintrf_mod, only : OBJ
use, non_intrinsic :: powalg_mod, only : hess_mul

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
real(RP), intent(inout) :: fval(:)  ! FVAL(NPT)
real(RP), intent(inout) :: gopt(:)  ! GOPT(N)
real(RP), intent(inout) :: hq(:, :)  ! HQ(N, N)
real(RP), intent(inout) :: pq(:)  ! PQ(NPT)
real(RP), intent(inout) :: sl(:)  ! SL(N)
real(RP), intent(inout) :: su(:)  ! SU(N)
real(RP), intent(inout) :: vlag(:)  ! VLAG(NPT+N)
real(RP), intent(inout) :: xbase(:)  ! XBASE(N)
real(RP), intent(inout) :: xhist(:, :)  ! XHIST(N, MAXXHIST)
real(RP), intent(inout) :: xopt(:)  ! XOPT(N)
real(RP), intent(inout) :: xpt(:, :)  ! XPT(N, NPT)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT-N-1)

! Outputs
!real(RP), intent(out) :: f
! Zaikun: 20220306
! For the moment, declare F as INOUT; otherwise, F is undefined if the subroutine exits in abnormal
! cases before the first evaluation of F. This should be resolved later when we can always keep the
! consistency between FOPT, XOPT, FVAL(KOPT), XPT(:, KOPT), F, and X. Maybe we do not even need to
! output F at all.
real(RP), intent(inout) :: f

! Local variables
character(len=*), parameter :: srname = 'RESCUE'
integer(IK) :: maxfhist
integer(IK) :: maxhist
integer(IK) :: maxxhist
integer(IK) :: n
integer(IK) :: npt
real(RP) :: ptsaux(2, size(xopt))
real(RP) :: ptsid(size(fval))
real(RP) :: x(size(xopt))
real(RP) :: beta, den(size(fval)), denom, moderr,      &
&        distsq(size(fval)), fbase, hdiag(size(fval)), sfrac,    &
&        temp, vlmxsq, vquad, dinc, xp, xq
integer(IK) :: ip, iq, iw, j, jp, jpn, k, &
&           kbase, knew, kold, kpt, np, nrem
real(RP) :: xpq(size(xopt)), pqinc(size(fval)), xxpt(size(fval)), wmv(size(xopt) + size(fval))
real(RP) :: bsum!, summ, vlag_test(size(xopt) + size(fval)), beta_test
logical :: mask(size(xopt))

n = int(size(xopt), kind(n))
npt = int(size(fval), kind(npt))
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
    call assert(size(vlag) == n + npt, 'SIZE(PQ) == N + NPT', srname)
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
    call assert(maxhist >= 0 .and. maxhist <= maxfun, '0 <= MAXHIST <= MAXFUN', srname)
end if


!
!     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
!       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
!       the corresponding arguments of BOBYQB on the entry to RESCUE.
!     NF is maintained as the number of calls of CALFUN so far, except that
!       NF is set to -1 if the value of MAXFUN prevents further progress.
!     KOPT is maintained so that FVAL(KOPT) is the least calculated function
!       value. Its correct value must be given on entry. It is updated if a
!       new least function value is found, but the corresponding changes to
!       XOPT and GOPT have to be made later by the calling program.
!     DELTA is the current trust region radius.
!     VLAG is a working space vector that will be used for the values of the
!       provisional Lagrange functions at each of the interpolation points.
!       They are part of a product that requires VLAG to be of length NDIM.
!     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
!       PTSAUX(2,J) specify the two positions of provisional interpolation
!       points when a nonzero step is taken along e_J (the J-th coordinate
!       direction) through XBASE+XOPT, as specified below. Usually these
!       steps have length DELTA, but other lengths are chosen if necessary
!       in order to satisfy the given bounds on the variables.
!     PTSID is also a working space array. It has NPT components that denote
!       provisional new positions of the original interpolation points, in
!       case changes are needed to restore the linear independence of the
!       interpolation conditions. The K-th point is a candidate for change
!       if and only if PTSID(K) is nonzero. In this case let p and q be the
!       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
!       and q are both positive, the step from XBASE+XOPT to the new K-th
!       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
!       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
!       p=0, respectively.
!     The first NDIM+NPT elements of the array W are used for working space.
!     The final elements of BMAT and ZMAT are set in a well-conditioned way
!       to the values that are appropriate for the new interpolation points.
!     The elements of GOPT, HQ and PQ are also revised to the values that are
!       appropriate to the final quadratic model.
!
!     Set some constants.
!

if (nf >= maxfun) then
    return
end if

kbase = kopt

np = n + 1
sfrac = HALF / real(np, RP)

! Shift the interpolation points so that XOPT becomes the origin, and set the elements of ZMAT to
! ZERO. The value of SUMPQ is required in the updating of HQ below. The squares of the distances
! from XOPT to the other interpolation points are set at DISTSQ. Increments of DINC may be added
! later to these squares to balance the consideration of the choice of point that is going to
! become current.
xpt = xpt - spread(xopt, dim=2, ncopies=npt)
distsq = sum((xpt)**2, dim=1)
dinc = maxval(distsq)
zmat = ZERO

! Update HQ so that HQ and PQ define the second derivatives of the model after XBASE has been
! shifted to the trust region centre.
xpq = matprod(xpt, pq) + HALF * sum(pq) * xopt
call r2update(hq, ONE, xopt, xpq)

! Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to ZERO, and set the elements of PTSAUX.
xbase = min(max(xl, xbase + xopt), xu)
sl = min(sl - xopt, ZERO)
su = max(su - xopt, ZERO)
xopt = ZERO
bmat = ZERO
ptsaux(1, :) = min(delta, su)
ptsaux(2, :) = max(-delta, sl)
mask = (ptsaux(1, :) + ptsaux(2, :) < 0)
ptsaux([1, 2], trueloc(mask)) = ptsaux([2, 1], trueloc(mask))
mask = (abs(ptsaux(2, :)) < HALF * abs(ptsaux(1, :)))
ptsaux(2, trueloc(mask)) = HALF * ptsaux(1, trueloc(mask))

fbase = fval(kopt)

! Set the identifiers of the artificial interpolation points that are along a coordinate direction
! from XOPT, and set the corresponding nonzero elements of BMAT and ZMAT.
ptsid(1) = sfrac
do j = 1, n
    jp = j + 1
    jpn = jp + n
    ptsid(jp) = real(j, RP) + sfrac
    if (jpn <= npt) then
        ptsid(jpn) = real(j, RP) / real(np, RP) + sfrac
        temp = ONE / (ptsaux(1, j) - ptsaux(2, j))
        bmat(j, jp) = -temp + ONE / ptsaux(1, j)
        bmat(j, jpn) = temp + ONE / ptsaux(2, j)
        bmat(j, 1) = -bmat(j, jp) - bmat(j, jpn)
        zmat(1, j) = sqrt(TWO) / abs(ptsaux(1, j) * ptsaux(2, j))
        zmat(jp, j) = zmat(1, j) * ptsaux(2, j) * temp
        zmat(jpn, j) = -zmat(1, j) * ptsaux(1, j) * temp
    else
        bmat(j, 1) = -ONE / ptsaux(1, j)
        bmat(j, jp) = ONE / ptsaux(1, j)
        bmat(j, j + npt) = -HALF * ptsaux(1, j)**2
    end if
end do

! Set any remaining identifiers with their nonzero elements of ZMAT.
if (npt >= n + np) then
    do k = 2 * np, npt
        iw = floor((real(k - np) - HALF) / real(n))
        ip = k - np - iw * n
        iq = ip + iw
        if (iq > n) iq = iq - n
        ptsid(k) = real(ip, RP) + real(iq, RP) / real(np, RP) + sfrac
        temp = ONE / (ptsaux(1, ip) * ptsaux(1, iq))
        zmat(1, k - np) = temp
        zmat(ip + 1, k - np) = -temp
        zmat(iq + 1, k - np) = -temp
        zmat(k, k - np) = temp
    end do
end if
nrem = npt
kold = 1
knew = kopt

80 continue

! Reorder the provisional points in the way that exchanges PTSID(KOLD) with PTSID(KNEW).
bmat(:, [kold, knew]) = bmat(:, [knew, kold])
zmat([kold, knew], :) = zmat([knew, kold], :)
ptsid(kold) = ptsid(knew)
ptsid(knew) = ZERO
distsq(knew) = ZERO
nrem = nrem - 1
if (knew /= kopt) then
    vlag([kold, knew]) = vlag([knew, kold])
    ! Update the BMAT and ZMAT matrices so that the status of the KNEW-th interpolation point can be
    ! changed from provisional to original. The branch to label 350 occurs if all the original
    ! points are reinstated. The nonnegative values of W(NDIM+K) are required in the search below.

    !----------------------------------------------------------------------------------------------!
    write (99, *) sum(zmat(knew, :)**2) * beta + vlag(knew)**2
    ! Without the writing, the following assertion may fail with NVFORTRAN.
    call assert(.not. abs(denom - (sum(zmat(knew, :)**2) * beta + vlag(knew)**2)) > 0, 'DENOM = DENOM_TEST', srname)
    !----------------------------------------------------------------------------------------------!
    call updateh(knew, beta, vlag, bmat, zmat)
    if (nrem == 0) goto 350
    distsq(1:npt) = abs(distsq(1:npt))
end if

120 continue

! Pick the index KNEW of an original interpolation point that has not yet replaced ONE of the
! provisional interpolation points, giving attention to the closeness to XOPT and to previous tries
! with KNEW.
if (any(distsq(1:npt) > 0)) then
    knew = int(minloc(distsq(1:npt), mask=(distsq(1:npt) > 0), dim=1), IK)
else
    goto 260
end if

! Calculate VLAG and BETA for the required updating of the H matrix if XPT(:, KNEW) is reinstated
! in the set of interpolation points.
! First, form the (W - V) vector for XPT(:, KNEW) as if it is going to replace XPT_TEST(:, KNEW)
! with the following XPT_TEST as defined in (5.4)--(5.5) of the BOBYQA paper for details.
! 1. XPT_TEST(:, KOPT) = 0;
! 2. For each K /= KOPT, if PTSID(K) == 0, then XPT_TEST(:, K) = XPT(:, K); if PTSID(K) > 0, then
! XPT_TEST(:, K) has nonzeros only at IP (if IP>0), IQ (if IQ>0) positions, where IP and IQ are
! the P(J) and Q(J) defined in and below (2.4) of the BOBYQA paper.
!
! In the code below, WMV = W - V = w(XNEW) - w(XOPT) without the (NPT+1)th entry, where
! XNEW = XPT(:,KNEW), XOPT = XPT_TEST = 0, and w(.) being defined by (6.3) of the NEWUOA paper (also
! well as (4.10) of the BOBYQA paper). Since XOPT = 0, we see from (6.3) that w(XOPT)(K) = 0 for all
! K except that w(XOPT)(NPT+1) = 1. Thus WMV is the same at w(XNEW) without the (NPT+1)-th entry.
! Therefore, WMV= [HALF*MATPROD(XNEW, XPT_TEST)**2, XNEW].
! Question (20220425): Can this be merged into CALVLAG/BETA in POWALG_MOD, maybe by modifying the
! signature of CALVLAG/BETA?
do k = 1, npt
    if (k == kopt) then
        wmv(k) = ZERO
    else if (ptsid(k) <= 0) then  ! Indeed, PTSID >= 0. So PTSID(K) <= 0 means PTSID(K) = 0.
        wmv(k) = inprod(xpt(:, knew), xpt(:, k))
    else
        ip = int(ptsid(k))  ! IP = 0 if 0 < PTSID(K) < 1. INT(X) = [X], i.e., it rounds X towards 0.
        iq = int(real(np, RP) * ptsid(k) - real(ip * np, RP))
        call assert(ip >= 0 .and. ip <= npt .and. iq >= 0 .and. iq <= npt, '0 <= IP, IQ <= NPT', srname)
        if (ip > 0 .and. iq > 0) then
            wmv(k) = xpt(ip, knew) * ptsaux(1, ip) + xpt(iq, knew) * ptsaux(1, iq)
        elseif (ip > 0) then
            wmv(k) = xpt(ip, knew) * ptsaux(1, ip)
        elseif (iq > 0) then
            wmv(k) = xpt(iq, knew) * ptsaux(2, iq)
        else
            wmv(k) = ZERO
        end if
    end if
    wmv(k) = HALF * wmv(k) * wmv(k)
end do
wmv(npt + 1:npt + n) = xpt(:, knew)
! Now calculate VLAG = H*WMV + e_KOPT according to (4.26) of the NEWUOA paper, except for VLAG(KOPT).
vlag(1:npt) = matprod(zmat, matprod(wmv(1:npt), zmat)) + matprod(wmv(npt + 1:npt + n), bmat(:, 1:npt))
vlag(npt + 1:npt + n) = matprod(bmat, wmv(1:npt + n))
! Now calculate BETA. According to (4.12) of the NEWUOA paper (as well as (4.10) of the BOBYQA paper),
! BETA = HALF*|XNEW - XOPT|^4 - WMV'*H*WMV. To calculate WMX'*H*WMV, note that
! WMV'*H*WMV = WMV' * [Z*Z', B2^T; B1, B2] * WMV with Z = ZMAT, B1 = BMAT(:, 1:NPT), and
! B2 = BMAT(:, NPT+1:NPT+N). Denoting W1 = WMV(1:NPT) and W2 = WMV(NPT+1:NPT+N), we then have
! WMV'*H*WMV = |W1*Z|^2 + 2*W1'*B1*W2 + W1'*B2*W2 = |W1*Z|^2 + W1'(B1*W2 + [B1, B2]*WMV).
bsum = inprod(wmv(1:n), matprod(bmat(:, 1:npt), wmv(1:npt)) + matprod(bmat, wmv))
beta = HALF * sum(xpt(:, knew)**2)**2 - sum(matprod(wmv(1:npt), zmat)**2) - bsum
! Finally, set VLAG(KOPT) to the correct value.
vlag(kopt) = vlag(kopt) + ONE

! KOLD is set to the index of the provisional interpolation point that is going to be deleted to
! make way for the KNEW-th original interpolation point. The choice of KOLD is governed by the
! avoidance of a small value of the denominator in the updating calculation of UPDATE.
hdiag = sum(zmat**2, dim=2)  ! Indeed, only HDIAG(PTSID /= 0) is needed.
!!MATLAB: hdiag(ptsid ~= 0) = sum(zmat(ptsid ~= 0, :), 2);
den = ZERO
where (ptsid /= 0) den = hdiag * beta + vlag(1:npt)**2
!!MATLAB: den(ptsid ~= 0) = hdiag(ptsid ~= 0)*beta + vlag(ptsid ~= 0)
kold = 0
denom = ZERO
if (any(den > 0)) then
    kold = int(maxloc(den, mask=(.not. is_nan(den)), dim=1), IK)
    denom = den(kold)
    !!MATLAB: [denom, kold] = max(den, [], 'omitnan');
end if
vlmxsq = maxval(vlag(1:npt)**2)
! KOLD > 0 is implied by DENOM > 1.0E-2*VLMXSQ (but NOT DENOM >= ...), yet we prefer to require
! KOLD > 0 explicitly.
if (kold > 0 .and. denom > 1.0E-2_RP * vlmxsq) then
    goto 80
else
    distsq(knew) = -distsq(knew) - dinc
    goto 120
end if

! All the final positions of the interpolation points have been chosen although any changes have not
! been included yet in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart from the
! shift of XBASE, the updating of the quadratic model remains to be done. The following cycle through
! the new interpolation points begins by putting the new point in XPT(:, KPT) and by setting PQ(KPT)
! to zero. A return occurs if MAXFUN prohibits another value of F or when all the new interpolation
! points are included in the model.
260 continue
do kpt = 1, npt
    if (ptsid(kpt) == ZERO) then
        cycle
    end if

    ! Absorb PQ(KPT)*XPT(:, KPT)*XPT(:, KPT)^T into the explicit part of the Hessian of the
    ! quadratic model. Implement R1UPDATE properly so that it ensures HQ is symmetric.
    call r1update(hq, pq(kpt), xpt(:, kpt))
    pq(kpt) = ZERO

    ip = int(ptsid(kpt))
    iq = int(real(np, RP) * ptsid(kpt) - real(ip * np, RP))

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
    ! multiply the KPT-th Lagrange function when the model is updated to provide interpolation to
    ! the new function value.
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
    if (is_nan(f) .or. is_posinf(f)) then
        exit
    end if
    if (f <= ftarget) then
        exit
    end if
    if (nf >= maxfun) then
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
    pq(trueloc(ptsid == 0)) = pq(trueloc(ptsid == 0)) + pqinc(trueloc(ptsid == 0))
    do k = 1, npt
        if (ptsid(k) == 0) then
            cycle
        end if
        ip = int(ptsid(k))
        iq = int(real(np, RP) * ptsid(k) - real(ip * np, RP))
        call assert(ip >= 0 .and. ip <= npt .and. iq >= 0 .and. iq <= npt, '0 <= IP, IQ <= NPT', srname)
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

if (kopt /= kbase) then
    xopt = xpt(:, kopt)
    gopt = gopt + hess_mul(xopt, xpt, pq, hq)
end if

350 continue

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
    call assert(size(vlag) == n + npt, 'SIZE(PQ) == N + NPT', srname)
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

return
end subroutine rescue


end module rescue_mod
