!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of bobyqb.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 11-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==bobyqb.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  SL,SU,XNEW,XALT,D,VLAG,W)
            subroutine BOBYQB(N, Npt, X, Xl, Xu, Rhobeg, Rhoend, Iprint,&
     & Maxfun, Xbase, Xpt, Fval, Xopt, Gopt, Hq, Pq, Bmat, Zmat, Ndim, S&
     &l, Su, Xnew, Xalt, D, Vlag, W, F, Info, Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--BOBYQB18
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            integer :: N
            integer :: Npt
            real*8, intent(INOUT), dimension(*) :: X
            real*8, dimension(*) :: Xl
            real*8, dimension(*) :: Xu
            real*8 :: Rhobeg
            real*8, intent(IN) :: Rhoend
            integer :: Iprint
            integer :: Maxfun
            real*8, intent(INOUT), dimension(*) :: Xbase
            real*8, intent(INOUT), dimension(Npt, *) :: Xpt
            real*8, intent(INOUT), dimension(*) :: Fval
            real*8, intent(INOUT), dimension(*) :: Xopt
            real*8, intent(INOUT), dimension(*) :: Gopt
            real*8, intent(INOUT), dimension(*) :: Hq
            real*8, intent(INOUT), dimension(*) :: Pq
            real*8, intent(INOUT), dimension(Ndim, *) :: Bmat
            real*8, dimension(Npt, *) :: Zmat
            integer :: Ndim
            real*8, intent(INOUT), dimension(*) :: Sl
            real*8, intent(INOUT), dimension(*) :: Su
            real*8, intent(INOUT), dimension(*) :: Xnew
            real*8, dimension(*) :: Xalt
            real*8, intent(INOUT), dimension(*) :: D
            real*8, intent(INOUT), dimension(*) :: Vlag
            real*8, intent(INOUT), dimension(*) :: W
            real*8, intent(INOUT) :: F
            integer, intent(OUT) :: Info
            real*8 :: Ftarget
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            real*8 :: adelt, almost_infinity, alpha, bdtest, bdtol, beta&
     &, biglsq, bsum, cauchy, crvmin, curv, delsq, delta, den, denom, de&
     &nsav, diff, diffa, diffb, diffc, dist, distsq, dnorm, dsq, dx, err&
     &big, fopt, fracsq, frhosq, gisq, gqsq, half, hdiag, one, pqold, ra&
     &tio, rho, scaden, sum, suma, sumb, sumpq, sumw, sumz, temp, ten, t&
     &enth, two, vquad, xoptsq, zero
            integer :: i, ih, ip, itest, j, jj, jp, k, kbase, knew, kopt&
     &, ksav, nf, nfsav, nh, np, nptm, nresc, ntrits
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
!       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT is a two-dimensional array that holds the coordinates of the
!       interpolation points relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XOPT is set to the displacement from XBASE of the trust region centre.
!     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
!       this factorization being ZMAT times ZMAT^T, which provides both the
!       correct rank and positive semi-definiteness.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
!       All the components of every XOPT are going to satisfy the bounds
!       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
!       XOPT is on a constraint boundary.
!     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
!       vector of variables for the next call of CALFUN. XNEW also satisfies
!       the SL and SU constraints in the way that has just been mentioned.
!     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
!       in order to increase the denominator in the updating of UPDATE.
!     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
!     VLAG contains the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     W is a one-dimensional array that is used for working space. Its length
!       must be at least 3*NDIM = 3*(NPT+N).
!
!     Set some constants.
!
            half = 0.5D0
            one = 1.0D0
            ten = 10.0D0
            tenth = 0.1D0
            two = 2.0D0
            zero = 0.0D0
            np = N + 1
            nptm = Npt - np
            nh = (N * np) / 2
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            almost_infinity = huge(0.0D0) / 2.D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, with the corresponding values of
!     of NF and KOPT, which are the number of calls of CALFUN so far and the
!     index of the interpolation point at the trust region centre. Then the
!     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
!     less than NPT. GOPT will be updated if KOPT is different from KBASE.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
            call PRELIM(N, Npt, X, Xl, Xu, Rhobeg, Iprint, Maxfun, Xbase&
     &, Xpt, Fval, Gopt, Hq, Pq, Bmat, Zmat, Ndim, Sl, Su, nf, kopt, F, &
     &Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            xoptsq = zero
            do i = 1, N
                Xopt(i) = Xpt(kopt, i)
                xoptsq = xoptsq + Xopt(i)**2
            end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: FSAVE is not needed any more. See line number 720.
!      FSAVE=FVAL(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom/Zaikun (on 04-06-2019/07-06-2019):
            if (F /= F .or. F > almost_infinity) then
                Info = -2
                goto 1000
            end if
!     By Tom (on 04-06-2019):
!     If F reached the target function, PRELIM will stop and BOBYQB
!     should stop here.
            if (F <= Ftarget) then
                Info = 1
                goto 1100
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (nf < Npt) then
                if (Iprint > 0) print 99007
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                Info = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                goto 1000
            end if
            kbase = 1
!
!     Complete the settings that are required for the iterative procedure.
!
            rho = Rhobeg
            delta = rho
            nresc = nf
            ntrits = 0
            diffa = zero
            diffb = zero
            itest = 0
            nfsav = nf
!
!     Update GOPT if necessary before the first iteration and after each
!     call of RESCUE that makes a call of CALFUN.
!
      100 if (kopt /= kbase) then
                ih = 0
                do j = 1, N
                    do i = 1, j
                        ih = ih + 1
                        if (i < j) Gopt(j) = Gopt(j) + Hq(ih) * Xopt(i)
                        Gopt(i) = Gopt(i) + Hq(ih) * Xopt(j)
                    end do
                end do
                if (nf > Npt) then
                    do k = 1, Npt
                        temp = zero
                        do j = 1, N
                            temp = temp + Xpt(k, j) * Xopt(j)
                        end do
                        temp = Pq(k) * temp
                        do i = 1, N
                            Gopt(i) = Gopt(i) + temp * Xpt(k, i)
                        end do
                    end do
                end if
            end if
!
!     Generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     The integer NTRITS is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. If the length
!     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
!     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRBOX, ALTMOV, or RESCUE is not predictable, and Segmentation Fault or
! infinite cycling may happen. This is because any equality/inequality
! comparison involving NaN returns FALSE, which can lead to unintended
! behavior of the code, including uninitialized indices.
!
!   60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,
      200 do i = 1, N
                if (Gopt(i) /= Gopt(i)) then
                    Info = -3
                    goto 1000
                end if
            end do
            do i = 1, nh
                if (Hq(i) /= Hq(i)) then
                    Info = -3
                    goto 1000
                end if
            end do
            do i = 1, Npt
                if (Pq(i) /= Pq(i)) then
                    Info = -3
                    goto 1000
                end if
            end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call TRSBOX(N, Npt, Xpt, Xopt, Gopt, Hq, Pq, Sl, Su, delta, &
     &Xnew, D, W, W(np), W(np + N), W(np + 2 * N), W(np + 3 * N), dsq, c&
     &rvmin)
            dnorm = DMIN1(delta, DSQRT(dsq))
            if (dnorm < half * rho) then
                ntrits = -1
                distsq = (ten * rho)**2
                if (nf <= nfsav + 2) goto 800
!
!     The following choice between labels 650 and 680 depends on whether or
!     not our work with the current RHO seems to be complete. Either RHO is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance HALF*RHO of XOPT.
!
                errbig = DMAX1(diffa, diffb, diffc)
                frhosq = 0.125D0 * rho * rho
                if (crvmin > zero .and. errbig > frhosq * crvmin) goto 8&
     &00
                bdtol = errbig / rho
                do j = 1, N
                    bdtest = bdtol
                    if (Xnew(j) == Sl(j)) bdtest = W(j)
                    if (Xnew(j) == Su(j)) bdtest = -W(j)
                    if (bdtest < bdtol) then
                        curv = Hq((j + j * j) / 2)
                        do k = 1, Npt
                            curv = curv + Pq(k) * Xpt(k, j)**2
                        end do
                        bdtest = bdtest + half * curv * rho
                        if (bdtest < bdtol) goto 800
                    end if
                end do
                goto 900
            end if
            ntrits = ntrits + 1
!
!     Severe cancellation is likely to occur if XOPT is too far from XBASE.
!     If the following test holds, then XBASE is shifted so that XOPT becomes
!     zero. The appropriate changes are made to BMAT and to the second
!     derivatives of the current model, beginning with the changes to BMAT
!     that do not depend on ZMAT. VLAG is used temporarily for working space.
!
      300 if (dsq <= 1.0D-3 * xoptsq) then
                fracsq = 0.25D0 * xoptsq
                sumpq = zero
                do k = 1, Npt
                    sumpq = sumpq + Pq(k)
                    sum = -half * xoptsq
                    do i = 1, N
                        sum = sum + Xpt(k, i) * Xopt(i)
                    end do
                    W(Npt + k) = sum
                    temp = fracsq - half * sum
                    do i = 1, N
                        W(i) = Bmat(k, i)
                        Vlag(i) = sum * Xpt(k, i) + temp * Xopt(i)
                        ip = Npt + i
                        do j = 1, i
                            Bmat(ip, j) = Bmat(ip, j) + W(i) * Vlag(j) +&
     & Vlag(i) * W(j)
                        end do
                    end do
                end do
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
                do jj = 1, nptm
                    sumz = zero
                    sumw = zero
                    do k = 1, Npt
                        sumz = sumz + Zmat(k, jj)
                        Vlag(k) = W(Npt + k) * Zmat(k, jj)
                        sumw = sumw + Vlag(k)
                    end do
                    do j = 1, N
                        sum = (fracsq * sumz - half * sumw) * Xopt(j)
                        do k = 1, Npt
                            sum = sum + Vlag(k) * Xpt(k, j)
                        end do
                        W(j) = sum
                        do k = 1, Npt
                            Bmat(k, j) = Bmat(k, j) + sum * Zmat(k, jj)
                        end do
                    end do
                    do i = 1, N
                        ip = i + Npt
                        temp = W(i)
                        do j = 1, i
                            Bmat(ip, j) = Bmat(ip, j) + temp * W(j)
                        end do
                    end do
                end do
!
!     The following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.
!
                ih = 0
                do j = 1, N
                    W(j) = -half * sumpq * Xopt(j)
                    do k = 1, Npt
                        W(j) = W(j) + Pq(k) * Xpt(k, j)
                        Xpt(k, j) = Xpt(k, j) - Xopt(j)
                    end do
                    do i = 1, j
                        ih = ih + 1
                        Hq(ih) = Hq(ih) + W(i) * Xopt(j) + Xopt(i) * W(j&
     &)
                        Bmat(Npt + i, j) = Bmat(Npt + j, i)
                    end do
                end do
                do i = 1, N
                    Xbase(i) = Xbase(i) + Xopt(i)
                    Xnew(i) = Xnew(i) - Xopt(i)
                    Sl(i) = Sl(i) - Xopt(i)
                    Su(i) = Su(i) - Xopt(i)
                    Xopt(i) = zero
                end do
                xoptsq = zero
            end if
            if (ntrits /= 0) goto 600
            goto 500
!
!     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
!     more expensive than the previous shift, because new matrices BMAT and
!     ZMAT are generated from scratch, which may include the replacement of
!     interpolation points whose positions seem to be causing near linear
!     dependence in the interpolation conditions. Therefore RESCUE is called
!     only if rounding errors have reduced by at least a factor of two the
!     denominator of the formula for updating the H matrix. It provides a
!     useful safeguard, but is not invoked in most applications of BOBYQA.
!
      400 nfsav = nf
            kbase = kopt
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29
! See the comments above line number 60.
            do i = 1, N
                if (Gopt(i) /= Gopt(i)) then
                    Info = -3
                    goto 1000
                end if
            end do
            do i = 1, nh
                if (Hq(i) /= Hq(i)) then
                    Info = -3
                    goto 1000
                end if
            end do
            do i = 1, Npt
                if (Pq(i) /= Pq(i)) then
                    Info = -3
                    goto 1000
                end if
            end do
            do j = 1, N
                do i = 1, Ndim
                    if (Bmat(i, j) /= Bmat(i, j)) then
                        Info = -3
                        goto 1000
                    end if
                end do
            end do
            do j = 1, nptm
                do i = 1, Npt
                    if (Zmat(i, j) /= Zmat(i, j)) then
                        Info = -3
                        goto 1000
                    end if
                end do
            end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  VLAG,W,W(N+NP),W(NDIM+NP))
            call RESCUE(N, Npt, Xl, Xu, Iprint, Maxfun, Xbase, Xpt, Fval&
     &, Xopt, Gopt, Hq, Pq, Bmat, Zmat, Ndim, Sl, Su, nf, delta, kopt, V&
     &lag, W, W(N + np), W(Ndim + np), F, Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     XOPT is updated now in case the branch below to label 720 is taken.
!     Any updating of GOPT occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.
!
            xoptsq = zero
            if (kopt /= kbase) then
                do i = 1, N
                    Xopt(i) = Xpt(kopt, i)
                    xoptsq = xoptsq + Xopt(i)**2
                end do
            end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom/Zaikun (on 04-06-2019/07-06-2019):
            if (F /= F .or. F > almost_infinity) then
                Info = -2
                goto 1000
            end if
!     By Tom (on 04-06-2019):
!     If F reached the target function, RESCUE will stop and BOBYQB
!     should stop here.
            if (F <= Ftarget) then
                Info = 1
                goto 1100
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (nf < 0) then
                nf = Maxfun
                if (Iprint > 0) print 99007
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                Info = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                goto 1000
            end if
            nresc = nf
            if (nfsav < nf) then
                nfsav = nf
                goto 100
            end if
            if (ntrits > 0) goto 200
!
!     Pick two alternative vectors of variables, relative to XBASE, that
!     are suitable as new positions of the KNEW-th interpolation point.
!     Firstly, XNEW is set to the point on a line through XOPT and another
!     interpolation point that minimizes the predicted value of the next
!     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
!     and SU bounds. Secondly, XALT is set to the best feasible point on
!     a constrained version of the Cauchy step of the KNEW-th Lagrange
!     function, the corresponding value of the square of this function
!     being returned in CAUCHY. The choice between these alternatives is
!     going to be made when the denominator is calculated.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  Zaikun 23-07-2019:
!  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
!     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
!
!  Although very rare, NaN can sometimes occur in BMAT or ZMAT. If it
!  happens, we terminate the code. See the comments above line number 60.
!  Indeed, if ALTMOV is called with such matrices, then altmov.f will
!  encounter a memory error at lines 173--174. This is because the first
!  value of PREDSQ in ALTOMOV (see line 159 of altmov.f) will be NaN, line
!  164 will not be reached, and hence no value will be assigned to IBDSAV.
!
!  Such an error was observed when BOBYQA was (mistakenly) tested on CUTEst
!  problem CONCON. CONCON is a nonlinearly constrained problem with
!  bounds. By mistake, BOBYQA was called to solve this problem,
!  neglecting all the constraints other than bounds. With only the bound
!  constraints, the objective function turned to be unbounded from
!  below, which led to abnormal values in BMAT (indeed, BETA defined in
!  lines 366--389 took NaN/infinite values).
!
      500 do j = 1, N
                do i = 1, Ndim
                    if (Bmat(i, j) /= Bmat(i, j)) then
                        Info = -3
                        goto 1000
                    end if
                end do
            end do
            do j = 1, nptm
                do i = 1, Npt
                    if (Zmat(i, j) /= Zmat(i, j)) then
                        Info = -3
                        goto 1000
                    end if
                end do
            end do
            call ALTMOV(N, Npt, Xpt, Xopt, Bmat, Zmat, Ndim, Sl, Su, kop&
     &t, knew, adelt, Xnew, Xalt, alpha, cauchy, W, W(np), W(Ndim + 1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do i = 1, N
                D(i) = Xnew(i) - Xopt(i)
            end do
      600 do
!
!     Calculate VLAG and BETA for the current choice of D. The scalar
!     product of D with XPT(K,.) is going to be held in W(NPT+K) for
!     use when VQUAD is calculated.
!
                do k = 1, Npt
                    suma = zero
                    sumb = zero
                    sum = zero
                    do j = 1, N
                        suma = suma + Xpt(k, j) * D(j)
                        sumb = sumb + Xpt(k, j) * Xopt(j)
                        sum = sum + Bmat(k, j) * D(j)
                    end do
                    W(k) = suma * (half * suma + sumb)
                    Vlag(k) = sum
                    W(Npt + k) = suma
                end do
                beta = zero
                do jj = 1, nptm
                    sum = zero
                    do k = 1, Npt
                        sum = sum + Zmat(k, jj) * W(k)
                    end do
                    beta = beta - sum * sum
                    do k = 1, Npt
                        Vlag(k) = Vlag(k) + sum * Zmat(k, jj)
                    end do
                end do
                dsq = zero
                bsum = zero
                dx = zero
                do j = 1, N
                    dsq = dsq + D(j)**2
                    sum = zero
                    do k = 1, Npt
                        sum = sum + W(k) * Bmat(k, j)
                    end do
                    bsum = bsum + sum * D(j)
                    jp = Npt + j
                    do i = 1, N
                        sum = sum + Bmat(jp, i) * D(i)
                    end do
                    Vlag(jp) = sum
                    bsum = bsum + sum * D(j)
                    dx = dx + D(j) * Xopt(j)
                end do
                beta = dx * dx + dsq * (xoptsq + dx + dx + half * dsq) +&
     & beta - bsum
                Vlag(kopt) = Vlag(kopt) + one
!
!     If NTRITS is zero, the denominator may be increased by replacing
!     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
!     rounding errors have damaged the chosen denominator.
!
                if (ntrits == 0) then
                    denom = Vlag(knew)**2 + alpha * beta
                    if (denom < cauchy .and. cauchy > zero) then
                        do i = 1, N
                            Xnew(i) = Xalt(i)
                            D(i) = Xnew(i) - Xopt(i)
                        end do
                        cauchy = zero
                        cycle
                    end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (DENOM .LE. HALF*VLAG(KNEW)**2) THEN
                    if (denom <= half * Vlag(knew)**2) then
!111111111111111111111!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if (nf > nresc) goto 400
                        if (Iprint > 0) print 99006
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                        Info = 4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        goto 1000
                    end if
!
!     Alternatively, if NTRITS is positive, then set KNEW to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. Again RESCUE may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     KNEW before calculating the next value of the objective function.
!
                else
                    delsq = delta * delta
                    scaden = zero
                    biglsq = zero
                    knew = 0
                    do k = 1, Npt
                        if (k == kopt) cycle
                        hdiag = zero
                        do jj = 1, nptm
                            hdiag = hdiag + Zmat(k, jj)**2
                        end do
                        den = beta * hdiag + Vlag(k)**2
                        distsq = zero
                        do j = 1, N
                            distsq = distsq + (Xpt(k, j) - Xopt(j))**2
                        end do
                        temp = DMAX1(one, (distsq / delsq)**2)
                        if (temp * den > scaden) then
                            scaden = temp * den
                            knew = k
                            denom = den
                        end if
                        biglsq = DMAX1(biglsq, temp * Vlag(k)**2)
                    end do
                    if (scaden <= half * biglsq) then
                        if (nf > nresc) goto 400
                        if (Iprint > 0) print 99006
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                        Info = 4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        goto 1000
                    end if
                end if
                exit
            end do
!
!     Put the variables for the next calculation of the objective function
!       in XNEW, with any adjustments for the bounds.
!
!
!     Calculate the value of the objective function at XBASE+XNEW, unless
!       the limit on the number of calculations of F has been reached.
!
      700 do i = 1, N
                X(i) = DMIN1(DMAX1(Xl(i), Xbase(i) + Xnew(i)), Xu(i))
                if (Xnew(i) == Sl(i)) X(i) = Xl(i)
                if (Xnew(i) == Su(i)) X(i) = Xu(i)
            end do
            if (nf >= Maxfun) then
                if (Iprint > 0) print 99007
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                Info = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                goto 1000
            end if
            nf = nf + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            do i = 1, N
                if (X(i) /= X(i)) then
                    F = X(i) ! Set F to NaN
                    if (nf == 1) then
                        fopt = F
                        Xopt(1:N) = zero
                    end if
                    Info = -1
                    goto 1000
                end if
            end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call CALFUN(N, X, F)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 04-06-2019):
            if (F /= F .or. F > almost_infinity) then
                if (nf == 1) then
                    fopt = F
                    Xopt(1:N) = zero
                end if
                Info = -2
                goto 1000
            end if
!     By Tom (on 04-06-2019):
!     If F achieves the function value, the algorithm exits.
            if (F <= Ftarget) then
                Info = 1
                goto 1100
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (Iprint == 3) then
                print 99001, nf, F, (X(i), i=1, N)
      99001 format(/4X, 'Function number', I6, ' F =', 1PD18.10, ' The c&
     &orresponding X is:'/(2X, 5D15.6))
            end if
            if (ntrits == -1) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: FSAVE is not needed any more. See line number 720.
!          FSAVE=F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                Info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                goto 1000
            end if
!
!     Use the quadratic model to predict the change in F due to the step D,
!       and set DIFF to the error of this prediction.
!
            fopt = Fval(kopt)
            vquad = zero
            ih = 0
            do j = 1, N
                vquad = vquad + D(j) * Gopt(j)
                do i = 1, j
                    ih = ih + 1
                    temp = D(i) * D(j)
                    if (i == j) temp = half * temp
                    vquad = vquad + Hq(ih) * temp
                end do
            end do
            do k = 1, Npt
                vquad = vquad + half * Pq(k) * W(Npt + k)**2
            end do
            diff = F - fopt - vquad
            diffc = diffb
            diffb = diffa
            diffa = DABS(diff)
            if (dnorm > rho) nfsav = nf
!
!     Pick the next value of DELTA after a trust region step.
!
            if (ntrits > 0) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (VQUAD .GE. ZERO) THEN
                if (vquad >= zero) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if (Iprint > 0) print 99002
      99002 format(/4X, 'Return from BOBYQA because a trust', ' region s&
     &tep has failed to reduce Q.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                    Info = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    goto 1000
                end if
                ratio = (F - fopt) / vquad
                if (ratio <= tenth) then
                    delta = DMIN1(half * delta, dnorm)
                elseif (ratio <= 0.7D0) then
                    delta = DMAX1(half * delta, dnorm)
                else
                    delta = DMAX1(half * delta, dnorm + dnorm)
                end if
                if (delta <= 1.5D0 * rho) delta = rho
!
!     Recalculate KNEW and DENOM if the new F is less than FOPT.
!
                if (F < fopt) then
                    ksav = knew
                    densav = denom
                    delsq = delta * delta
                    scaden = zero
                    biglsq = zero
                    knew = 0
                    do k = 1, Npt
                        hdiag = zero
                        do jj = 1, nptm
                            hdiag = hdiag + Zmat(k, jj)**2
                        end do
                        den = beta * hdiag + Vlag(k)**2
                        distsq = zero
                        do j = 1, N
                            distsq = distsq + (Xpt(k, j) - Xnew(j))**2
                        end do
                        temp = DMAX1(one, (distsq / delsq)**2)
                        if (temp * den > scaden) then
                            scaden = temp * den
                            knew = k
                            denom = den
                        end if
                        biglsq = DMAX1(biglsq, temp * Vlag(k)**2)
                    end do
                    if (scaden <= half * biglsq) then
                        knew = ksav
                        denom = densav
                    end if
                end if
            end if
!
!     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
!     moved. Also update the second derivative terms of the model.
!
            call UPDATE(N, Npt, Bmat, Zmat, Ndim, Vlag, beta, denom, kne&
     &w, W)
            ih = 0
            pqold = Pq(knew)
            Pq(knew) = zero
            do i = 1, N
                temp = pqold * Xpt(knew, i)
                do j = 1, i
                    ih = ih + 1
                    Hq(ih) = Hq(ih) + temp * Xpt(knew, j)
                end do
            end do
            do jj = 1, nptm
                temp = diff * Zmat(knew, jj)
                do k = 1, Npt
                    Pq(k) = Pq(k) + temp * Zmat(k, jj)
                end do
            end do
!
!     Include the new interpolation point, and make the changes to GOPT at
!     the old XOPT that are caused by the updating of the quadratic model.
!
            Fval(knew) = F
            do i = 1, N
                Xpt(knew, i) = Xnew(i)
                W(i) = Bmat(knew, i)
            end do
            do k = 1, Npt
                suma = zero
                do jj = 1, nptm
                    suma = suma + Zmat(knew, jj) * Zmat(k, jj)
                end do
                sumb = zero
                do j = 1, N
                    sumb = sumb + Xpt(k, j) * Xopt(j)
                end do
                temp = suma * sumb
                do i = 1, N
                    W(i) = W(i) + temp * Xpt(k, i)
                end do
            end do
            do i = 1, N
                Gopt(i) = Gopt(i) + diff * W(i)
            end do
!
!     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
!
            if (F < fopt) then
                kopt = knew
                xoptsq = zero
                ih = 0
                do j = 1, N
                    Xopt(j) = Xnew(j)
                    xoptsq = xoptsq + Xopt(j)**2
                    do i = 1, j
                        ih = ih + 1
                        if (i < j) Gopt(j) = Gopt(j) + Hq(ih) * D(i)
                        Gopt(i) = Gopt(i) + Hq(ih) * D(j)
                    end do
                end do
                do k = 1, Npt
                    temp = zero
                    do j = 1, N
                        temp = temp + Xpt(k, j) * D(j)
                    end do
                    temp = Pq(k) * temp
                    do i = 1, N
                        Gopt(i) = Gopt(i) + temp * Xpt(k, i)
                    end do
                end do
            end if
!
!     Calculate the parameters of the least Frobenius norm interpolant to
!     the current data, the gradient of this interpolant at XOPT being put
!     into VLAG(NPT+I), I=1,2,...,N.
!
            if (ntrits > 0) then
                do k = 1, Npt
                    Vlag(k) = Fval(k) - Fval(kopt)
                    W(k) = zero
                end do
                do j = 1, nptm
                    sum = zero
                    do k = 1, Npt
                        sum = sum + Zmat(k, j) * Vlag(k)
                    end do
                    do k = 1, Npt
                        W(k) = W(k) + sum * Zmat(k, j)
                    end do
                end do
                do k = 1, Npt
                    sum = zero
                    do j = 1, N
                        sum = sum + Xpt(k, j) * Xopt(j)
                    end do
                    W(k + Npt) = W(k)
                    W(k) = sum * W(k)
                end do
                gqsq = zero
                gisq = zero
                do i = 1, N
                    sum = zero
                    do k = 1, Npt
                        sum = sum + Bmat(k, i) * Vlag(k) + Xpt(k, i) * W&
     &(k)
                    end do
                    if (Xopt(i) == Sl(i)) then
                        gqsq = gqsq + DMIN1(zero, Gopt(i))**2
                        gisq = gisq + DMIN1(zero, sum)**2
                    elseif (Xopt(i) == Su(i)) then
                        gqsq = gqsq + DMAX1(zero, Gopt(i))**2
                        gisq = gisq + DMAX1(zero, sum)**2
                    else
                        gqsq = gqsq + Gopt(i)**2
                        gisq = gisq + sum * sum
                    end if
                    Vlag(Npt + i) = sum
                end do
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
                itest = itest + 1
                if (gqsq < ten * gisq) itest = 0
                if (itest >= 3) then
                    do i = 1, MAX0(Npt, nh)
                        if (i <= N) Gopt(i) = Vlag(Npt + i)
                        if (i <= Npt) Pq(i) = W(Npt + i)
                        if (i <= nh) Hq(i) = zero
                        itest = 0
                    end do
                end if
            end if
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case NTRITS=0 occurs
!     when the new interpolation point was reached by an alternative step.
!
            if (ntrits == 0) goto 200
            if (F <= fopt + tenth * vquad) goto 200
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
            distsq = DMAX1((two * delta)**2, (ten * rho)**2)
      800 knew = 0
            do k = 1, Npt
                sum = zero
                do j = 1, N
                    sum = sum + (Xpt(k, j) - Xopt(j))**2
                end do
                if (sum > distsq) then
                    knew = k
                    distsq = sum
                end if
            end do
!
!     If KNEW is positive, then ALTMOV finds alternative new positions for
!     the KNEW-th interpolation point within distance ADELT of XOPT. It is
!     reached via label 90. Otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current RHO are complete.
!
            if (knew > 0) then
                dist = DSQRT(distsq)
                if (ntrits == -1) then
                    delta = DMIN1(tenth * delta, half * dist)
                    if (delta <= 1.5D0 * rho) delta = rho
                end if
                ntrits = 0
                adelt = DMAX1(DMIN1(tenth * dist, delta), rho)
                dsq = adelt * adelt
                goto 300
            end if
            if (ntrits /= -1) then
                if (ratio > zero) goto 200
                if (DMAX1(delta, dnorm) > rho) goto 200
            end if
!
!     The calculations with the current value of RHO are complete. Pick the
!       next values of RHO and DELTA.
!
      900 if (rho > Rhoend) then
                delta = half * rho
                ratio = rho / Rhoend
                if (ratio <= 16.0D0) then
                    rho = Rhoend
                elseif (ratio <= 250.0D0) then
                    rho = DSQRT(ratio) * Rhoend
                else
                    rho = tenth * rho
                end if
                delta = DMAX1(delta, rho)
                if (Iprint >= 2) then
                    if (Iprint >= 3) print 99003
      99003 format(5X)
                    print 99004, rho, nf
      99004 format(/4X, 'New RHO =', 1PD11.4, 5X, 'Number of', ' functio&
     &n values =', I6)
                    print 99008, Fval(kopt), (Xbase(i) + Xopt(i), i=1, N&
     &)
                end if
                ntrits = 0
                nfsav = nf
                goto 200
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            else
                Info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
!
!     Return from the calculation, after another Newton-Raphson step, if
!       it is too short to have been tried before.
!
            if (ntrits == -1) goto 700
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  720 IF (FVAL(KOPT) .LE. FSAVE) THEN
!  Why update X only when FVAL(KOPT) .LE. FSAVE? This seems INCORRECT,
!  because it may lead to a return with F and X that are not the best
!  available.
      1000 if (Fval(kopt) <= F .or. F /= F) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                do i = 1, N
                    X(i) = DMIN1(DMAX1(Xl(i), Xbase(i) + Xopt(i)), Xu(i)&
     &)
                    if (Xopt(i) == Sl(i)) X(i) = Xl(i)
                    if (Xopt(i) == Su(i)) X(i) = Xu(i)
                end do
                F = Fval(kopt)
            end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (IPRINT .GE. 1) THEN
      1100 if (Iprint >= 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                print 99005, nf
      99005 format(/4X, 'At the return from BOBYQA', 5X, 'Number of func&
     &tion values =', I6)
                print 99008, F, (X(i), i=1, N)
            end if
      99006 format(/5X, 'Return from BOBYQA because of much', ' cancella&
     &tion in a denominator.')
      99007 format(/4X, 'Return from BOBYQA because CALFUN has been', ' &
     &called MAXFUN times.')
      99008 format(4X, 'Least value of F =', 1PD23.15, 9X, 'The correspo&
     &nding X is:'/(2X, 5D15.6))
            end subroutine BOBYQB