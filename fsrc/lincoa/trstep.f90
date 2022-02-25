subroutine trstep(n, npt, m, amat, xpt, hq, pq, nact, iact, rescon, &
     &  qfac, rfac, snorm, step, g, resnew, resact, d, dw, w)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK

! Solver-specific modules
use, non_intrinsic :: getact_mod, only : getact

implicit real(RP) (a - h, o - z)
implicit integer(IK) (i - n)
dimension amat(n, *), xpt(npt, *), hq(*), pq(*), iact(*), &
& rescon(*), qfac(n, *), rfac(*), step(*), g(*), resnew(*), resact(*), &
& d(*), dw(*), w(*)
!
!     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
!       are the same as the terms with these names in LINCOB. If RESCON(J)
!       is negative, then |RESCON(J)| must be no less than the trust region
!       radius, so that the J-th constraint can be ignored.
!     SNORM is set to the trust region radius DELTA initially. On the
!       return, however, it is the length of the calculated STEP, which is
!       set to zero if the constraints do not allow a long enough step.
!     STEP is the total calculated step so far from the trust region centre,
!       its final value being given by the sequence of CG iterations, which
!       terminate if the trust region boundary is reached.
!     G must be set on entry to the gradient of the quadratic model at the
!       trust region centre. It is used as working space, however, and is
!       always the gradient of the model at the current STEP, except that
!       on return the value of G(1) is set to ONE instead of to ZERO if
!       and only if GETACT is called more than once.
!     RESNEW, RESACT, D, DW and W are used for working space. A negative
!       value of RESNEW(J) indicates that the J-th constraint does not
!       restrict the CG steps of the current trust region calculation, a
!       zero value of RESNEW(J) indicates that the J-th constraint is active,
!       and otherwise RESNEW(J) is set to the greater of TINY and the actual
!       residual of the J-th constraint for the current STEP. RESACT holds
!       the residuals of the active constraints, which may be positive.
!       D is the search direction of each line search. DW is either another
!       search direction or the change in gradient along D. The length of W
!       must be at least MAX[M,2*N].
!
!     Set some numbers for the conjugate gradient iterations.
!
half = 0.5D0
one = 1.0D0
tinynum = real(tiny(0.0), RP)
zero = 0.0D0
ctest = 0.01D0
snsq = snorm * snorm
!
!     Set the initial elements of RESNEW, RESACT and STEP.
!
if (m > 0) then
    do j = 1, m
        resnew(j) = rescon(j)
        if (rescon(j) >= snorm) then
            resnew(j) = -one
        else if (rescon(j) >= zero) then
            resnew(j) = max(resnew(j), tinynum)
        end if
    end do
    if (nact > 0) then
        do k = 1, nact
            resact(k) = rescon(iact(k))
            resnew(iact(k)) = zero
        end do
    end if
end if
do i = 1, n
    step(i) = zero
end do
ss = zero
reduct = zero
ncall = 0
!
!     GETACT picks the active set for the current STEP. It also sets DW to
!       the vector closest to -G that is orthogonal to the normals of the
!       active constraints. DW is scaled to have length 0.2*SNORM, as then
!       a move of DW from STEP is allowed by the linear constraints.
!
40 ncall = ncall + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: B is never used in GETACT
!      CALL GETACT (N,M,AMAT,B,NACT,IACT,QFAC,RFAC,SNORM,RESNEW,
call getact(n, m, amat, nact, iact, qfac, rfac, snorm, resnew, &
& resact, g, dw, w, w(n + 1))
if (w(n + 1) == zero) goto 320
scale = 0.2D0 * snorm / sqrt(w(n + 1))
do i = 1, n
    dw(i) = scale * dw(i)
end do
!
!     If the modulus of the residual of an active constraint is substantial,
!       then set D to the shortest move from STEP to the boundaries of the
!       active constraints.
!
resmax = zero
if (nact > 0) then
    do k = 1, nact
        resmax = max(resmax, resact(k))
    end do
end if
gamma = zero
if (resmax > 1.0D-4 * snorm) then
    ir = 0
    do k = 1, nact
        temp = resact(k)
        if (k >= 2) then
            do i = 1, k - 1
                ir = ir + 1
                temp = temp - rfac(ir) * w(i)
            end do
        end if
        ir = ir + 1
        w(k) = temp / rfac(ir)
    end do
    do i = 1, n
        d(i) = zero
        do k = 1, nact
            d(i) = d(i) + w(k) * qfac(i, k)
        end do
    end do
!
!     The vector D that has just been calculated is also the shortest move
!       from STEP+DW to the boundaries of the active constraints. Set GAMMA
!       to the greatest steplength of this move that satisfies the trust
!       region bound.
!
    rhs = snsq
    ds = zero
    dd = zero
    do i = 1, n
        sum = step(i) + dw(i)
        rhs = rhs - sum * sum
        ds = ds + d(i) * sum
        dd = dd + d(i)**2
    end do
    if (rhs > zero) then
        temp = sqrt(ds * ds + dd * rhs)
        if (ds <= zero) then
            ! Zaikun 20210925
            ! What if we are at the first iteration? BLEN = DELTA/|D|? See TRSAPP.F90 of NEWUOA.
            gamma = (temp - ds) / dd
        else
            gamma = rhs / (temp + ds)
        end if
    end if
!
!     Reduce the steplength GAMMA if necessary so that the move along D
!       also satisfies the linear constraints.
!
    j = 0
110 if (gamma > zero) then
        j = j + 1
        if (resnew(j) > zero) then
            ad = zero
            adw = zero
            do i = 1, n
                ad = ad + amat(i, j) * d(i)
                adw = adw + amat(i, j) * dw(i)
            end do
            if (ad > zero) then
                temp = max((resnew(j) - adw) / ad, zero)
                gamma = min(gamma, temp)
            end if
        end if
        if (j < m) goto 110
    end if
    gamma = min(gamma, one)
end if
!
!     Set the next direction for seeking a reduction in the model function
!       subject to the trust region bound and the linear constraints.
!
if (gamma <= zero) then
    do i = 1, n
        d(i) = dw(i)
    end do
    icount = nact
else
    do i = 1, n
        d(i) = dw(i) + gamma * d(i)
    end do
    icount = nact - 1
end if
alpbd = one
!
!     Set ALPHA to the steplength from STEP along D to the trust region
!       boundary. Return if the first derivative term of this step is
!       sufficiently small or if no further progress is possible.
!
150 icount = icount + 1
rhs = snsq - ss
if (rhs <= zero) goto 320
dg = zero
ds = zero
dd = zero
do i = 1, n
    dg = dg + d(i) * g(i)
    ds = ds + d(i) * step(i)
    dd = dd + d(i)**2
end do
if (dg >= zero) goto 320
temp = sqrt(rhs * dd + ds * ds)
if (ds <= zero) then
    alpha = (temp - ds) / dd
else
    alpha = rhs / (temp + ds)
end if
if (-alpha * dg <= ctest * reduct) goto 320
!
!     Set DW to the change in gradient along D.
!
ih = 0
do j = 1, n
    dw(j) = zero
    do i = 1, j
        ih = ih + 1
        if (i < j) dw(j) = dw(j) + hq(ih) * d(i)
        dw(i) = dw(i) + hq(ih) * d(j)
    end do
end do
do k = 1, npt
    temp = zero
    do j = 1, n
        temp = temp + xpt(k, j) * d(j)
    end do
    temp = pq(k) * temp
    do i = 1, n
        dw(i) = dw(i) + temp * xpt(k, i)
    end do
end do
!
!     Set DGD to the curvature of the model along D. Then reduce ALPHA if
!       necessary to the value that minimizes the model.
!
dgd = zero
do i = 1, n
    dgd = dgd + d(i) * dw(i)
end do
alpht = alpha
if (dg + alpha * dgd > zero) then
    alpha = -dg / dgd
end if
!
!     Make a further reduction in ALPHA if necessary to preserve feasibility,
!       and put some scalar products of D with constraint gradients in W.
!
alphm = alpha
jsav = 0
if (m > 0) then
    do j = 1, m
        ad = zero
        if (resnew(j) > zero) then
            do i = 1, n
                ad = ad + amat(i, j) * d(i)
            end do
            if (alpha * ad > resnew(j)) then
                alpha = resnew(j) / ad
                jsav = j
            end if
        end if
        w(j) = ad
    end do
end if
alpha = max(alpha, alpbd)
alpha = min(alpha, alphm)
if (icount == nact) alpha = min(alpha, one)
!
!     Update STEP, G, RESNEW, RESACT and REDUCT.
!
ss = zero
do i = 1, n
    step(i) = step(i) + alpha * d(i)
    ss = ss + step(i)**2
    g(i) = g(i) + alpha * dw(i)
end do
if (m > 0) then
    do j = 1, m
        if (resnew(j) > zero) then
            resnew(j) = max(resnew(j) - alpha * w(j), tinynum)
        end if
    end do
end if
if (icount == nact .and. nact > 0) then
    do k = 1, nact
        resact(k) = (one - gamma) * resact(k)
    end do
end if
reduct = reduct - alpha * (dg + half * alpha * dgd)
!
!     Test for termination. Branch to label 40 if there is a new active
!       constraint and if the distance from STEP to the trust region
!       boundary is at least 0.2*SNORM.
!
! Zaikun 2019-08-29: the code can encounter infinite cycling due to NaN
! values. Exit when NCALL is large or NaN detected.
if (ncall > min(10000, 100 * (m + 1) * n) .or.  &
& alpha /= alpha .or. alpht /= alpht .or. &
& alphm /= alphm .or. dgd /= dgd .or. dg /= dg .or. &
& ss /= ss .or. snsq /= snsq .or. reduct /= reduct) then
    goto 320
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (alpha == alpht) goto 320
temp = -alphm * (dg + half * alphm * dgd)
if (temp <= ctest * reduct) goto 320
if (jsav > 0) then
    if (ss <= 0.64D0 * snsq) goto 40
    goto 320
end if
if (icount == n) goto 320
!
!     Calculate the next search direction, which is conjugate to the
!       previous one except in the case ICOUNT=NACT.
!
if (nact > 0) then
    do j = nact + 1, n
        w(j) = zero
        do i = 1, n
            w(j) = w(j) + g(i) * qfac(i, j)
        end do
    end do
    do i = 1, n
        temp = zero
        do j = nact + 1, n
            temp = temp + qfac(i, j) * w(j)
        end do
        w(n + i) = temp
    end do
else
    do i = 1, n
        w(n + i) = g(i)
    end do
end if
if (icount == nact) then
    beta = zero
else
    wgd = zero
    do i = 1, n
        wgd = wgd + w(n + i) * dw(i)
    end do
    beta = wgd / dgd
end if
do i = 1, n
    d(i) = -w(n + i) + beta * d(i)
end do
alpbd = zero
goto 150
!
!     Return from the subroutine.
!
320 snorm = zero
if (reduct > zero) snorm = sqrt(ss)
g(1) = zero
if (ncall > 1) g(1) = one
return
end
