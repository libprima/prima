!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of rescue.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 06-Jul-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==rescue.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  KOPT,VLAG,PTSAUX,PTSID,W)
            subroutine RESCUE(N, Npt, Xl, Xu, Iprint, Maxfun, Xbase, Xpt&
     &, Fval, Xopt, Gopt, Hq, Pq, Bmat, Zmat, Ndim, Sl, Su, Nf, Delta, K&
     &opt, Vlag, Ptsaux, Ptsid, W, F, Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--RESCUE14
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            integer :: N
            integer :: Npt
            real*8, intent(IN), dimension(*) :: Xl
            real*8, intent(IN), dimension(*) :: Xu
            integer, intent(IN) :: Iprint
            integer, intent(IN) :: Maxfun
            real*8, intent(INOUT), dimension(*) :: Xbase
            real*8, intent(INOUT), dimension(Npt, *) :: Xpt
            real*8, intent(INOUT), dimension(*) :: Fval
            real*8, intent(INOUT), dimension(*) :: Xopt
            real*8, intent(INOUT), dimension(*) :: Gopt
            real*8, intent(INOUT), dimension(*) :: Hq
            real*8, intent(INOUT), dimension(*) :: Pq
            real*8, intent(INOUT), dimension(Ndim, *) :: Bmat
            real*8, intent(INOUT), dimension(Npt, *) :: Zmat
            integer :: Ndim
            real*8, intent(INOUT), dimension(*) :: Sl
            real*8, intent(INOUT), dimension(*) :: Su
            integer, intent(INOUT) :: Nf
            real*8, intent(IN) :: Delta
            integer, intent(INOUT) :: Kopt
            real*8, intent(INOUT), dimension(*) :: Vlag
            real*8, intent(INOUT), dimension(2, *) :: Ptsaux
            real*8, intent(INOUT), dimension(*) :: Ptsid
            real*8, intent(INOUT), dimension(*) :: W
            real*8 :: F
            real*8, intent(IN) :: Ftarget
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            real*8 :: almost_infinity, beta, bsum, den, denom, diff, dis&
     &tsq, dsqmin, fbase, half, hdiag, one, sfrac, sum, sumpq, temp, vlm&
     &xsq, vquad, winc, xp, xq, zero
            integer :: i, ih, ihp, ihq, ip, iq, iw, j, jp, jpn, k, knew,&
     & kold, kpt, np, nptm, nrem
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            half = 0.5D0
            one = 1.0D0
            zero = 0.0D0
            np = N + 1
            sfrac = half / DFLOAT(np)
            nptm = Npt - np
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            almost_infinity = huge(0.0D0) / 2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Shift the interpolation points so that XOPT becomes the origin, and set
!     the elements of ZMAT to zero. The value of SUMPQ is required in the
!     updating of HQ below. The squares of the distances from XOPT to the
!     other interpolation points are set at the end of W. Increments of WINC
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.
!
            sumpq = zero
            winc = zero
            do k = 1, Npt
                distsq = zero
                do j = 1, N
                    Xpt(k, j) = Xpt(k, j) - Xopt(j)
                    distsq = distsq + Xpt(k, j)**2
                end do
                sumpq = sumpq + Pq(k)
                W(Ndim + k) = distsq
                winc = DMAX1(winc, distsq)
                do j = 1, nptm
                    Zmat(k, j) = zero
                end do
            end do
!
!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.
!
            ih = 0
            do j = 1, N
                W(j) = half * sumpq * Xopt(j)
                do k = 1, Npt
                    W(j) = W(j) + Pq(k) * Xpt(k, j)
                end do
                do i = 1, j
                    ih = ih + 1
                    Hq(ih) = Hq(ih) + W(i) * Xopt(j) + W(j) * Xopt(i)
                end do
            end do
!
!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
!     also set the elements of PTSAUX.
!
            do j = 1, N
                Xbase(j) = Xbase(j) + Xopt(j)
                Sl(j) = Sl(j) - Xopt(j)
                Su(j) = Su(j) - Xopt(j)
                Xopt(j) = zero
                Ptsaux(1, j) = DMIN1(Delta, Su(j))
                Ptsaux(2, j) = DMAX1(-Delta, Sl(j))
                if (Ptsaux(1, j) + Ptsaux(2, j) < zero) then
                    temp = Ptsaux(1, j)
                    Ptsaux(1, j) = Ptsaux(2, j)
                    Ptsaux(2, j) = temp
                end if
                if (DABS(Ptsaux(2, j)) < half * DABS(Ptsaux(1, j))) Ptsa&
     &ux(2, j) = half * Ptsaux(1, j)
                do i = 1, Ndim
                    Bmat(i, j) = zero
                end do
            end do
            fbase = Fval(Kopt)
!
!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonzero elements of BMAT and ZMAT.
!
            Ptsid(1) = sfrac
            do j = 1, N
                jp = j + 1
                jpn = jp + N
                Ptsid(jp) = DFLOAT(j) + sfrac
                if (jpn <= Npt) then
                    Ptsid(jpn) = DFLOAT(j) / DFLOAT(np) + sfrac
                    temp = one / (Ptsaux(1, j) - Ptsaux(2, j))
                    Bmat(jp, j) = -temp + one / Ptsaux(1, j)
                    Bmat(jpn, j) = temp + one / Ptsaux(2, j)
                    Bmat(1, j) = -Bmat(jp, j) - Bmat(jpn, j)
                    Zmat(1, j) = DSQRT(2.0D0) / DABS(Ptsaux(1, j) * Ptsa&
     &ux(2, j))
                    Zmat(jp, j) = Zmat(1, j) * Ptsaux(2, j) * temp
                    Zmat(jpn, j) = -Zmat(1, j) * Ptsaux(1, j) * temp
                else
                    Bmat(1, j) = -one / Ptsaux(1, j)
                    Bmat(jp, j) = one / Ptsaux(1, j)
                    Bmat(j + Npt, j) = -half * Ptsaux(1, j)**2
                end if
            end do
!
!     Set any remaining identifiers with their nonzero elements of ZMAT.
!
            if (Npt >= N + np) then
                do k = 2 * np, Npt
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
                    iw = int((DFLOAT(k - np) - half) / DFLOAT(N))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ip = k - np - iw * N
                    iq = ip + iw
                    if (iq > N) iq = iq - N
                    Ptsid(k) = DFLOAT(ip) + DFLOAT(iq) / DFLOAT(np) + sf&
     &rac
                    temp = one / (Ptsaux(1, ip) * Ptsaux(1, iq))
                    Zmat(1, k - np) = temp
                    Zmat(ip + 1, k - np) = -temp
                    Zmat(iq + 1, k - np) = -temp
                    Zmat(k, k - np) = temp
                end do
            end if
            nrem = Npt
            kold = 1
            knew = Kopt
!
!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).
!
100   do j = 1, N
                temp = Bmat(kold, j)
                Bmat(kold, j) = Bmat(knew, j)
                Bmat(knew, j) = temp
            end do
            do j = 1, nptm
                temp = Zmat(kold, j)
                Zmat(kold, j) = Zmat(knew, j)
                Zmat(knew, j) = temp
            end do
            Ptsid(kold) = Ptsid(knew)
            Ptsid(knew) = zero
            W(Ndim + knew) = zero
            nrem = nrem - 1
            if (knew /= Kopt) then
                temp = Vlag(kold)
                Vlag(kold) = Vlag(knew)
                Vlag(knew) = temp
!
!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.
!
                call UPDATE(N, Npt, Bmat, Zmat, Ndim, Vlag, beta, denom,&
     & knew, W)
                if (nrem == 0) goto 99999
                do k = 1, Npt
                    W(Ndim + k) = DABS(W(Ndim + k))
                end do
            end if
            do
!
!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.
!
                dsqmin = zero
                do k = 1, Npt
                    if (W(Ndim + k) > zero) then
                        if (dsqmin == zero .or. W(Ndim + k) < dsqmin) th&
     &en
                            knew = k
                            dsqmin = W(Ndim + k)
                        end if
                    end if
                end do
                if (dsqmin == zero) then
!
!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be done. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.
!
                    do kpt = 1, Npt
                        if (Ptsid(kpt) == zero) cycle
                        if (Nf >= Maxfun) then
                            Nf = -1
                            exit
                        end if
                        ih = 0
                        do j = 1, N
                            W(j) = Xpt(kpt, j)
                            Xpt(kpt, j) = zero
                            temp = Pq(kpt) * W(j)
                            do i = 1, j
                                ih = ih + 1
                                Hq(ih) = Hq(ih) + temp * W(i)
                            end do
                        end do
                        Pq(kpt) = zero
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IP=PTSID(KPT)
!      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
                        ip = int(Ptsid(kpt))
                        iq = int(DFLOAT(np) * Ptsid(kpt) - DFLOAT(ip * n&
     &p))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if (ip > 0) then
                            xp = Ptsaux(1, ip)
                            Xpt(kpt, ip) = xp
                        end if
                        if (iq > 0) then
                            xq = Ptsaux(1, iq)
                            if (ip == 0) xq = Ptsaux(2, iq)
                            Xpt(kpt, iq) = xq
                        end if
!
!     Set VQUAD to the value of the current model at the new point.
!
                        vquad = fbase
                        if (ip > 0) then
                            ihp = (ip + ip * ip) / 2
                            vquad = vquad + xp * (Gopt(ip) + half * xp *&
     & Hq(ihp))
                        end if
                        if (iq > 0) then
                            ihq = (iq + iq * iq) / 2
                            vquad = vquad + xq * (Gopt(iq) + half * xq *&
     & Hq(ihq))
                            if (ip > 0) then
                                iw = MAX0(ihp, ihq) - IABS(ip - iq)
                                vquad = vquad + xp * xq * Hq(iw)
                            end if
                        end if
                        do k = 1, Npt
                            temp = zero
                            if (ip > 0) temp = temp + xp * Xpt(k, ip)
                            if (iq > 0) temp = temp + xq * Xpt(k, iq)
                            vquad = vquad + half * Pq(k) * temp * temp
                        end do
!
!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
                        do i = 1, N
                            W(i) = DMIN1(DMAX1(Xl(i), Xbase(i) + Xpt(kpt&
     &, i)), Xu(i))
                            if (Xpt(kpt, i) == Sl(i)) W(i) = Xl(i)
                            if (Xpt(kpt, i) == Su(i)) W(i) = Xu(i)
                        end do
                        Nf = Nf + 1
                        call CALFUN(N, W, F)
                        if (Iprint == 3) then
                            print 99001, Nf, F, (W(i), i=1, N)
99001 format(/4X, 'Function number', I6, ' F =', 1PD18.10, ' The corresp&
     &onding X is:'/(2X, 5D15.6))
                        end if
                        Fval(kpt) = F
                        if (F < Fval(Kopt)) Kopt = kpt
                        diff = F - vquad
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
                        do i = 1, N
                            Gopt(i) = Gopt(i) + diff * Bmat(kpt, i)
                        end do
                        do k = 1, Npt
                            sum = zero
                            do j = 1, nptm
                                sum = sum + Zmat(k, j) * Zmat(kpt, j)
                            end do
                            temp = diff * sum
                            if (Ptsid(k) == zero) then
                                Pq(k) = Pq(k) + temp
                            else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IP=PTSID(K)
!          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
                                ip = int(Ptsid(k))
                                iq = int(DFLOAT(np) * Ptsid(k) - DFLOAT(&
     &ip * np))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                ihq = (iq * iq + iq) / 2
                                if (ip == 0) then
                                    Hq(ihq) = Hq(ihq) + temp * Ptsaux(2,&
     & iq)**2
                                else
                                    ihp = (ip * ip + ip) / 2
                                    Hq(ihp) = Hq(ihp) + temp * Ptsaux(1,&
     & ip)**2
                                    if (iq > 0) then
                                        Hq(ihq) = Hq(ihq) + temp * Ptsau&
     &x(1, iq)**2
                                        iw = MAX0(ihp, ihq) - IABS(iq - &
     &ip)
                                        Hq(iw) = Hq(iw) + temp * Ptsaux(&
     &1, ip) * Ptsaux(1, iq)
                                    end if
                                end if
                            end if
                        end do
                        Ptsid(kpt) = zero
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 03-06-2019):
!     If a NaN or an infinite value has been reached during the
!     evaluation of the objective function, the loop exit after setting
!     all the parameters, not to raise an exception. KOPT is set to KPT
!     to check in BOBYQB weather FVAL(KOPT) is NaN or infinite value or
!     not.
                        if (F /= F .or. F > almost_infinity) exit
!     By Tom (on 04-06-2019):
!     If the target function value is reached, the loop exit and KOPT is
!     set to KPT to check in BOBYQB weather FVAL(KOPT) .LE. FTARGET
                        if (F <= Ftarget) exit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end do
                    exit
                else
!
!     Form the W-vector of the chosen original interpolation point.
!
                    do j = 1, N
                        W(Npt + j) = Xpt(knew, j)
                    end do
                    do k = 1, Npt
                        sum = zero
                        if (k == Kopt) then
                        elseif (Ptsid(k) == zero) then
                            do j = 1, N
                                sum = sum + W(Npt + j) * Xpt(k, j)
                            end do
                        else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IP=PTSID(K)
                            ip = int(Ptsid(k))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            if (ip > 0) sum = W(Npt + ip) * Ptsaux(1, ip&
     &)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
                            iq = int(DFLOAT(np) * Ptsid(k) - DFLOAT(ip *&
     & np))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            if (iq > 0) then
                                iw = 1
                                if (ip == 0) iw = 2
                                sum = sum + W(Npt + iq) * Ptsaux(iw, iq)
                            end if
                        end if
                        W(k) = half * sum * sum
                    end do
!
!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.
!
                    do k = 1, Npt
                        sum = zero
                        do j = 1, N
                            sum = sum + Bmat(k, j) * W(Npt + j)
                        end do
                        Vlag(k) = sum
                    end do
                    beta = zero
                    do j = 1, nptm
                        sum = zero
                        do k = 1, Npt
                            sum = sum + Zmat(k, j) * W(k)
                        end do
                        beta = beta - sum * sum
                        do k = 1, Npt
                            Vlag(k) = Vlag(k) + sum * Zmat(k, j)
                        end do
                    end do
                    bsum = zero
                    distsq = zero
                    do j = 1, N
                        sum = zero
                        do k = 1, Npt
                            sum = sum + Bmat(k, j) * W(k)
                        end do
                        jp = j + Npt
                        bsum = bsum + sum * W(jp)
                        do ip = Npt + 1, Ndim
                            sum = sum + Bmat(ip, j) * W(ip)
                        end do
                        bsum = bsum + sum * W(jp)
                        Vlag(jp) = sum
                        distsq = distsq + Xpt(knew, j)**2
                    end do
                    beta = half * distsq * distsq + beta - bsum
                    Vlag(Kopt) = Vlag(Kopt) + one
!
!     KOLD is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the KNEW-th original interpolation
!     point. The choice of KOLD is governed by the avoidance of a small value
!     of the denominator in the updating calculation of UPDATE.
!
                    denom = zero
                    vlmxsq = zero
                    do k = 1, Npt
                        if (Ptsid(k) /= zero) then
                            hdiag = zero
                            do j = 1, nptm
                                hdiag = hdiag + Zmat(k, j)**2
                            end do
                            den = beta * hdiag + Vlag(k)**2
                            if (den > denom) then
                                kold = k
                                denom = den
                            end if
                        end if
                        vlmxsq = DMAX1(vlmxsq, Vlag(k)**2)
                    end do
                    if (denom <= 1.0D-2 * vlmxsq) then
                        W(Ndim + knew) = -W(Ndim + knew) - winc
                        cycle
                    end if
                    goto 100
                end if
            end do
99999 end subroutine RESCUE