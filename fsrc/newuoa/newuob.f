      subroutine newuob (n, npt, x, rhobeg, rhoend, iprint, maxfun,     &
     & xbase, xopt, xnew, xpt, fval, gq, hq, pq, bmat, zmat, ndim,      &
     & d, vlag, w, f, info, ftarget)

      use pdfomod, only : rp, zero, one, half, tenth, is_nan, is_posinf
      implicit none

      ! inputs
      integer, intent(in) :: n, npt, iprint, maxfun, ndim
      integer, intent(out) :: info
      real(kind = rp), intent(in) :: rhobeg, rhoend, ftarget
      real(kind = rp), intent(out) :: f
      real(kind = rp), intent(inout) :: x(n), xbase(n), xopt(n),        &
     & xnew(n), xpt(npt, n), fval(npt), gq(n), hq((n*(n+1))/2), pq(npt)
      real(kind = rp), intent(inout) :: bmat(npt + n, n),               &
     & zmat(npt, npt - n - 1), d(n), vlag(npt + n), w(10*(npt + n))

      ! other variables
      integer :: i, idz, ih, itest, j, k, knew, kopt, ksav, ktemp, nf,  &
     & nfsav, nftest, subinfo
      real(kind = rp) :: alpha, beta, crvmin, delta, detrat, diff,      &
     & diffa, diffb, diffc, distsq, dnorm, dsq, dstep
      real(kind = rp) :: fopt, fsave, gisq, gqsq, hdiag, ratio, rho,    &
     & rhosq, summation, temp, vquad, xoptsq, wcheck(npt + n)

      ! The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
      ! identical to the corresponding arguments in SUBROUTINE NEWUOA.
      ! XBASE will hold a shift of origin that should reduce the
      ! contributions from rounding errors to values of the model and
      ! Lagrange functions.
      ! XOPT will be set to the displacement from XBASE of the vector of
      ! variables that provides the least calculated F so far.
      ! XNEW will be set to the displacement from XBASE of the vector of
      ! variables for the current calculation of F.
      ! XPT will contain the interpolation point coordinates relative to
      ! XBASE.
      ! FVAL will hold the values of F at the interpolation points.
      ! GQ will hold the gradient of the quadratic model at XBASE.
      ! HQ will hold the explicit second order derivatives of the
      ! quadratic model.
      ! PQ will contain the parameters of the implicit second order
      ! derivatives of the quadratic model.
      ! BMAT will hold the last N columns of H. ZMAT will hold the
      ! factorization of the leading NPT by NPT submatrix of H, this
      ! factorization being ZMAT times Diag(DZ) times ZMAT^T, where the
      ! elements of DZ are plus or minus one, as specified by IDZ.
      ! NDIM is the first dimension of BMAT and has the value NPT + N.
      ! D is reserved for trial steps from XOPT.
      ! VLAG will contain the values of the Lagrange functions at a new
      ! point X. They are part of a product that requires VLAG to be of
      ! length NDIM = NPT+N.
      ! The array W will be used for working space. Its length must be
      ! at least 10*NDIM = 10*(NPT + N).

      ! Set some constants.
      nftest = max(maxfun, 1)

      call initialize(n, npt, rhobeg, x, xbase, xpt, f, fval, xopt,     &
     & fopt, kopt, bmat, zmat, gq, hq, pq, nf, subinfo, ftarget)
      if (subinfo == 1 .or. subinfo == -1 .or. subinfo == -2 .or.       &
     & subinfo == -3) then
          info = subinfo
          goto 530
      end if

      ! Set some more initial values.
      rho = rhobeg
      delta = rho
      idz = 1
      diffa = zero
      diffb = zero
      itest = 0
      nfsav = nf
      xopt = xpt(kopt, :)
      xoptsq = zero
      do i = 1, n
          xoptsq = xoptsq + xopt(i)**2
      end do

      ! Begin the iterative procedure.

      ! Generate the next trust region step and test its length. Set
      ! KNEW to -1 if the purpose of the next F is to improve the model.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur
      ! in the models. In such a case, we terminate the code. Otherwise,
      ! the behavior of TRSAPP, BIGDEN, or BIGLAG is not predictable,
      ! and Segmentation Fault or infinite cycling may happen. This is
      ! because any equality/inequality comparison involving NaN returns
      ! FALSE, which can lead to unintended behavior of the code,
      ! including uninitialized indices.
  100 if (any(is_nan(gq)) .or. any(is_nan(hq)) .or. any(is_nan(pq)))then
          info = -3
          goto 530
      end if

      knew = 0

      call trsapp (n, npt, xopt, xpt, gq, hq, pq, delta, d, w, w(n+1),  &
     & w(2*n+1), w(3*n+1), crvmin)
      !dsq = dot_product(d, d)
      dsq = zero
      do i = 1, n
          dsq = dsq + d(i)**2
      end do
      dnorm = dmin1(delta, sqrt(dsq))
      if (dnorm < half*rho) then
          knew = -1
          delta = tenth*delta
          ratio = -one
          if (delta <= 1.5_rp*rho) delta = rho
          if (0.125_rp*crvmin*rho*rho <= max(diffa, diffb, diffc) .or.  &
     &     nf <= nfsav + 2) then 
              goto 460
          else
              goto 490
          end if
      end if

      ! Shift XBASE if XOPT may be too far from XBASE.
      if (dsq <= 1.0e-3_rp*xoptsq) then
          call shiftbase(n, npt, idz, xopt, pq, bmat, zmat, gq, hq, xpt)
          xbase = xbase + xopt
          xopt = zero
          xoptsq = zero
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Zaikun 2019-08-29: See the comments below line number 100
  120 if (any(is_nan(bmat)) .or. any(is_nan(zmat)) .or. any(is_nan(gq)) &
     & .or. any(is_nan(hq)) .or. any(is_nan(pq))) then
          info = -3
          goto 530
      end if
!  120 do j = 1, n
!          do i = 1, ndim
!              if (bmat(i, j) /= bmat(i, j)) then
!                  info = -3
!                  goto 530
!              end if
!          end do
!      end do
!      do j = 1, npt - n - 1
!          do i = 1, npt
!              if (zmat(i, j) /= zmat(i, j)) then
!                  info = -3
!                  goto 530
!              end if
!          end do
!      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Pick the model step if KNEW is positive. A different choice of
      ! D may be made later, if the choice of D by BIGLAG causes
      ! substantial cancellation in DENOM.
      if (knew > 0) then
          call biglag (n, npt, xopt, xpt, bmat, zmat, idz, ndim, knew,  &
     &     dstep, d, alpha, vlag, vlag(npt + 1), w, w(n+1), w(2*n+1))
      end if

      ! Calculate VLAG and BETA for the current choice of D. The first
      ! NPT components of W_check will be held in W.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if ((dsq-sum(d**2))/dsq > 1e-15) print *, (dsq-sum(d**2))/dsq 
!      if(knew > 0 .and. dsq /= dstep*dstep) print *, knew, dsq-dstep**2 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! NOTE:
      ! 1. w(1:npt) (i.e., wcheck) may be used later 
      ! 2. dsq is supposed to be ||d||^2, but it is calculated by
      !    dstep^2 when knew > 0. This is exact in theory, because dstep
      !    is the trust region radius in biglag, where ||d|| should equal
      !    this radius at return. However, in computation, there may be
      !    a slight difference between the two due to rounding errors.
      !    On the other hand, xoptsq is always equal to ||xopt||^2. 

      !dsq = dot_product(d, d)
!      call vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d, vlag,  &
!     & beta, wcheck)
      call vlagbeta(n, npt, idz, kopt, bmat, zmat, xpt, xopt, d, vlag,  &
     & beta, wcheck, dsq, xoptsq)
      !w(1:npt) = wcheck  

      ! If KNEW is positive and if the cancellation in DENOM is
      ! unacceptable, then BIGDEN calculates an alternative model step,
      ! XNEW being used for working space.
      if (knew > 0) then
          if (abs(one + alpha*beta/vlag(knew)**2) <= 0.8_rp) then
              call bigden (n, npt, xopt, xpt, bmat, zmat, idz, ndim,    &
     &         kopt, knew, d, wcheck, vlag, beta, xnew,w(ndim+1),       &
     &         w(6*ndim+1))
          end if
      end if
      
      ! Calculate the next value of the objective function.
      xnew = xopt + d
      x = xbase + xnew
      if (any(is_nan(x))) then
          f = sum(x)  ! Set F to NaN. It is necessary.
          info = -1
          goto 530
      else
          call calfun(n, x, f)
          nf = nf + 1
      end if

      if (iprint == 3) then
          print 330, nf, f, (x(i), i = 1, n)
  330     FORMAT (/4X, 'Function number', I6, '    F = ', 1PD18.10,     &
     &     '    The corresponding X is:'/(2X, 5D15.6))
      end if
      
      ! Use the quadratic model to predict the change in F due to the
      ! step D,  and set DIFF to the error of this prediction.
      vquad = zero
      ih = 0
      do j = 1, n
          vquad = vquad + d(j)*gq(j)
          do i = 1, j
              ih = ih + 1
              temp = d(i)*xnew(j) + d(j)*xopt(i)
              if (i == j) temp = half*temp
              vquad = vquad + temp*hq(ih)
          end do
      end do
      do k = 1, npt
          vquad = vquad + pq(k)*wcheck(k)
      end do
      diff = f - fopt - vquad

      diffc = diffb
      diffb = diffa
      diffa = abs(diff)

      if (dnorm > rho) nfsav = nf

      ! Update FOPT and XOPT if the new F is the least value of the
      ! objective function so far. The branch when KNEW is positive
      ! occurs if D is not a trust region step.
      fsave = fopt
      if (f < fopt) then
          fopt = f
          xopt = xnew
          xoptsq = zero
          do i = 1, n
              xoptsq = xoptsq + xopt(i)**2
          end do
      end if
      ksav = knew
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! By Zaikun (commented on 02-06-2019; implemented in 2016):
      ! Return if F has an NaN or almost infinite value.
      ! If this happends at the first function evaluation (i.e., NF=1),
      ! then it is necessary to set FOPT and XOPT before going to 530,
      ! because these two variables have not been set yet (line 70
      ! will not be reached).
      if (is_nan(f) .or. is_posinf(f)) then 
          info = -2
          goto 530
      end if
      ! By Zaikun (commented on 02-06-2019; implemented in 2016):
      ! Exit if F .LE. FTARGET.
      if (f <= ftarget) then
          info = 1
          goto 530
      end if
      if (nf >= nftest) then
          if (iprint > 0) print 320
  320         FORMAT (/4X, 'Return from NEWUOA because CALFUN has       &
     & been called MAXFUN times.')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          info = 3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          goto 530
      end if

!      if (knew == -1) goto 530  !!!

      if (knew > 0) goto 410

      ! Pick the next value of DELTA after a trust region step.
      ! IF (VQUAD .GE. ZERO) THEN
      if (is_nan(vquad) .or. vquad >= zero) then
          if (iprint > 0) print 370
  370         FORMAT (/4X, 'Return from NEWUOA because a trust          &
     & region step has failed to reduce Q.')
          info = 2
          goto 530
      end if
      ratio = (f - fsave)/vquad
      if (ratio <= tenth) then
          delta = half*dnorm
      else if (ratio <= 0.7_rp) then
          delta = max(half*delta, dnorm)
      else
          delta = max(half*delta, dnorm + dnorm)
      end if
      if (delta <= 1.5_rp*rho) delta = rho
      
      ! Set KNEW to the index of the next interpolation point to delete.
      rhosq = max(tenth*delta, rho)**2
      ktemp = 0
      detrat = zero
      if (f >= fsave) then
          ktemp = kopt
          detrat = one
      end if
      do k = 1, npt
          hdiag = zero
          do j = 1, npt - n - 1
              temp = one
              if (j < idz) temp = -one
              hdiag = hdiag + temp*zmat(k, j)**2
          end do
          temp = abs(beta*hdiag + vlag(k)**2)
          distsq = zero
          do j = 1, n
              distsq = distsq + (xpt(k, j) - xopt(j))**2
          end do
          if (distsq > rhosq) temp = temp*(distsq/rhosq)**3
          if (temp > detrat .and. k /= ktemp) then
              detrat = temp
              knew = k
          end if
      end do
      if (knew == 0) goto 460
      !
      ! Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation
      ! point can be moved. Begin the updating of the quadratic model,
      ! starting with the explicit second derivative term.
      !
  410 call update (n, npt, bmat, zmat, idz, ndim, vlag, beta,knew,w)
      fval(knew) = f
      ih = 0
      do i = 1, n
          temp = pq(knew)*xpt(knew, i)
          do j = 1, i
              ih = ih + 1
              hq(ih) = hq(ih) + temp*xpt(knew, j)
          end do
      end do
      pq(knew) = zero
      !
      ! Update the other second derivative parameters, and then the
      ! gradient of the model. Also include the new interpolation point.
      do j = 1, npt - n - 1
          temp = diff*zmat(knew, j)
          if (j < idz) temp = -temp
          do k = 1, npt
              pq(k) = pq(k) + temp*zmat(k, j)
          end do
      end do
      xpt(knew, :) = xnew
      gq = gq + diff*bmat(knew, :)
      gqsq = zero
      do i = 1, n
          gqsq = gqsq + gq(i)**2
      end do

      ! If a trust region step makes a small change to the objective
      ! function, then calculate the gradient of the least Frobenius
      ! norm interpolant at XBASE, and store it in W, using VLAG for
      ! a vector of right hand sides.
      if (ksav == 0 .and. delta == rho) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Zaikun 2019-08-26: It is observed in Zhang Zaikun's PhD
          ! thesis (Section 3.3.2) that it is more reasonable and more
          ! efficient to check the value of RATIO instead of ABS(RATIO).
          ! IF (DABS(RATIO) .GT. 1.0D-2) THEN
          if (ratio > 1.0e-2_rp) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              itest = 0
          else
              do k = 1, npt
                  vlag(k) = fval(k) - fval(kopt)
              end do
              gisq = zero
              do i = 1, n
                  summation = zero
                  do k = 1, npt
                      summation = summation + bmat(k, i)*vlag(k)
                  end do
                  gisq = gisq + summation*summation
                  w(i) = summation
              end do

              ! Test whether to replace the new quadratic model by the
              ! least Frobenius norm interpolant, making the replacement
              ! if the test is satisfied.
              itest = itest + 1
              if (gqsq < 100.0_rp*gisq) itest = 0
              if (itest >= 3) then
                  gq = w(1:n)
                  hq = zero
                  do j = 1, npt - n - 1
                      w(j) = zero
                      do k = 1, npt
                          w(j) = w(j) + vlag(k)*zmat(k, j)
                      end do
                      if (j < idz) w(j) = -w(j)
                  end do
                  pq = zero
                  do k = 1, npt
                      do j = 1, npt - n - 1
                          pq(k) = pq(k) + zmat(k, j)*w(j)
                      end do
                  end do
                  itest = 0
              end if
          end if
      end if
      if (f < fsave) kopt = knew
    
      ! If a trust region step has provided a sufficient decrease in F,
      ! then branch for another trust region calculation. The case
      ! KSAVE>0 occurs when the new function value was calculated by
      ! a model step.
      if (f <= fsave + tenth*vquad) goto 100
      if (ksav > 0) goto 100
      
      ! Alternatively, find out if the interpolation points are close
      ! enough to the best point so far.
      knew = 0

  460 distsq = 4.0_rp*delta*delta
      do k = 1, npt
          summation = zero
          do j = 1, n
              summation = summation + (xpt(k, j) - xopt(j))**2
          end do
          if (summation > distsq) then
              knew = k
              distsq = summation
          end if
      end do

      ! If KNEW is positive, then set DSTEP, and branch back for the
      ! next iteration, which will generate a "model step".
      if (knew > 0) then
          dstep = max(dmin1(tenth*sqrt(distsq), half*delta), rho)
          dsq = dstep*dstep
          if (dsq <= 1.0e-3_rp*xoptsq) then
              call shiftbase(n, npt, idz, xopt, pq, bmat, zmat, gq, hq, &
     &         xpt)
              xbase = xbase + xopt
              xopt = zero
              xoptsq = zero
          end if
          goto 120
      end if

      if (ratio > zero .or. max(delta, dnorm) > rho) goto 100
      
      ! The calculations with the current value of RHO are complete.
      ! Pick the next values of RHO and DELTA.
  490 if (rho > rhoend) then
          delta = half*rho
          ratio = rho/rhoend
          if (ratio <= 16.0_rp) then
              rho = rhoend
          else if (ratio <= 250.0_rp) then
              rho = sqrt(ratio)*rhoend
          else
              rho = tenth*rho
          end if
          delta = max(delta, rho)
          if (iprint >= 2) then
              if (iprint >= 3) print 500
  500             FORMAT (5X)
              print 510, rho, nf
  510         FORMAT (/4X, 'New RHO = ', 1PD11.4, 5X,                   &
     &         'Number of function values = ', I6)
              print 520, fopt, (xbase + xopt)
  520         FORMAT (4X, 'Least value of F = ', 1PD23.15, 9X,          &
     &         'The corresponding X is:'/(2X, 5D15.6))
          end if

          nfsav = nf

          goto 100 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
          info = 0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
      
      ! Return from the calculation, after another Newton-Raphson step,
      ! if it is too short to have been tried before.
      if (knew == -1) then
          x = xbase + (xopt + d)
          if (any(is_nan(x))) then
              f = sum(x)  ! Set F to NaN. It is necessary.
              info = -1
          else
              call calfun(n, x, f)
              nf = nf + 1
          end if
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! By Zaikun (commented on 02-06-2019; implemented in 2016):
      ! Note that (FOPT .LE. F) is FALSE if F is NaN; When F is NaN, it
      ! is also necessary to update X and F.
      ! 530 IF (FOPT .LE. F) THEN
  530 if (is_nan(f) .or. fopt <= f) then
          x = xbase + xopt
          f = fopt
      end if
      if (iprint >= 1) then
          print 550, nf
  550     FORMAT (/4X, 'At the return from NEWUOA', 5X,                 &
     &     'Number of function values = ', I6)
          print 520, f, (x(i), i = 1, n)
      end if
      return

      end subroutine newuob


CCCCCCCCCCCCCCCCCCCCCCCCCCC Auxillary Subroutines CCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE QALT(GQ, HQ, PQ, FVAL, SMAT, ZMAT, N, NPT, NPTM, KOPT, 
     +    IDZ)
C
C     QALT calculates the alternative model, namely the model that
C     minimizes the F-norm of the Hessian subject to the interpolation
C     conditions. 
C
C     Note that SMAT = BMAT(1:NPT, 1:N)
C
          IMPLICIT NONE
          INTEGER, PARAMETER :: DP = KIND(0.0D0) 
          ! DP IS THE KIND FOR DOUBLE PRECISION
          INTEGER, INTENT(IN) :: N, NPT, NPTM, KOPT, IDZ
          REAL(KIND = DP), INTENT(IN) :: FVAL(NPT), SMAT(NPT, N), 
     +    ZMAT(NPT, NPTM)
          REAL(KIND = DP), INTENT(OUT) :: GQ(N), HQ(N*(N+1)/2), PQ(NPT)
          REAL(KIND = DP) :: VLAG(NPT), W(NPTM)

          VLAG = FVAL - FVAL(KOPT)
          GQ = MATMUL(VLAG, SMAT)
          HQ = 0.0D0
          W = MATMUL(VLAG, ZMAT)
          W(1:IDZ-1) = - W(1:IDZ-1)
          PQ = MATMUL(ZMAT, W)

          RETURN

      END SUBROUTINE QALT


      SUBROUTINE DQ(QDIFF, D, X, XPT, GQ, HQ, PQ, N, NPT)
C
C     DQ calculates QDIFF = Q(Y) - Q(X), where Y = X + D
C
          IMPLICIT NONE
          INTEGER, PARAMETER :: DP = KIND(0.0D0) 
          ! DP IS THE KIND FOR DOUBLE PRECISION
          INTEGER, INTENT(IN) :: N, NPT
          REAL(KIND = DP), INTENT(IN) :: D(N), X(N), XPT(NPT, N), 
     +    GQ(N), HQ(N*(N+1)/2), PQ(NPT)
          REAL(KIND = DP), INTENT(OUT) :: QDIFF
          INTEGER :: I, J, IH
          REAL(KIND = DP) :: S(N), SD

          ! S = X + Y 
          S = D + X + X 

          ! 1st order term plus explicit 2nd-order term
          QDIFF = DOT_PRODUCT(D, GQ) + 
     +            0.5D0*SUM(PQ*MATMUL(XPT,S)*MATMUL(XPT,D)) 

          ! Implicit 2nd-order term
          IH = 0 
          DO I = 1, N
              DO J = 1, I 
                  IF (I .EQ. J) THEN
                      SD = S(I)*D(I)
                  ELSE
                      SD = S(I)*D(J) + S(J)*D(I)
                  END IF
                  IH = IH + 1
                  QDIFF = QDIFF + 0.5D0 * HQ(IH) * SD 
              END DO
          END DO

          RETURN

      END SUBROUTINE DQ


C      SUBROUTINE TESTINT(EINT, FVAL, XPT, GQ, HQ, PQ, N, NPT, KOPT)
CC
CC     TESTINT tests how well Q interpolates F.
CC     At return, EINT(K) is the interpolation error at XPT(K, :)
CC
C          IMPLICIT NONE
C          INTEGER, PARAMETER :: DP = KIND(0.0D0) 
C          ! DP IS THE KIND FOR DOUBLE PRECISION
C          INTEGER, INTENT(IN) :: N, NPT, KOPT
C          REAL(KIND = DP), INTENT(IN) :: FVAL(NPT), XPT(NPT, N), GQ(N), 
C     +    HQ(N*(N+1)/2), PQ(NPT)
C          REAL(KIND = DP), INTENT(OUT) :: EINT(NPT) 
C          INTEGER :: K
C          REAL(KIND = DP) :: FOPT, FREF, XOPT(N), D(N), QDIFF, FDIFF
C
C          FOPT = FVAL(KOPT)
C          FREF = MAX(1.0D0, MAXVAL(ABS(FVAL-FOPT)))
C          XOPT = XPT(KOPT, :)
C          EINT(KOPT) = 0.0D0
C          DO K = 1, NPT
C              IF (K .NE. KOPT) THEN
C                  D = XPT(K, :) - XOPT
C                  CALL DQ(QDIFF, D, XOPT, XPT, GQ, HQ, PQ, N, NPT)
C                  FDIFF = FVAL(K) - FOPT
C                  EINT(K) = ABS(FDIFF - QDIFF)/FREF
C              END IF
C          END DO
C        
C          RETURN
C
C      END SUBROUTINE
