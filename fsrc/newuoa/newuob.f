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
      integer :: i, idz, itest, k, knew, kopt, ksav, nf,                &
     & nfsav, nftest, subinfo
      real(kind = rp) :: alpha, beta, crvmin, delta, prederr(3),        &
     & distsq, dnorm, dsq, dstep, xdiff(n), xdsq(npt)
      real(kind = rp) :: fopt, fsave, galt(n), galtsq, gqsq, hdiag(npt),&
     & ratio, rho, rhosq, vquad, xoptsq, wcheck(npt+n), sigma(npt) 

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
      prederr = zero
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
          if (0.125_rp*crvmin*rho*rho <= maxval(abs(prederr)) .or.      &
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
      ! step D,  and set PREDERR to the error of this prediction.
      
      !call calquad(vquad, d, xopt, xpt, gq, hq, pq, n, npt)
      call calquad(vquad, d, xopt, xpt, gq, hq, pq, n,npt,wcheck(1:npt))

      prederr(2 : size(prederr)) = prederr(1 : size(prederr)-1)
      prederr(1) = f - fopt - vquad

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
      do k = 1, npt
          hdiag(k) = -sum(zmat(k, 1 : idz - 1)**2) +                    &
     &     sum(zmat(k, idz : npt - n - 1)**2)
          xdiff = xpt(k, :) - xopt
          xdsq(k) = dot_product(xdiff, xdiff) 
      end do
      sigma = abs(beta*hdiag + vlag(1 : npt)**2)
      sigma = sigma * max(xdsq/rhosq, one)**3
      if (f >= fsave) then
      ! Set SIGMA(KOPT) = -1 to prevent KNEW from being KOPT 
          sigma(kopt) = -one  
      end if 
      if (maxval(sigma) > one .or. f < fsave) then
      ! KNEW > 0 unless maxval(sigma) <= 1 and f >= fsave
          knew = maxloc(sigma, 1)
      else
          knew = 0
      end if

      if (knew == 0) goto 460
      
      ! Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation
      ! point can be moved. Begin the updating of the quadratic model,
      ! starting with the explicit second derivative term.
  410 call update (n, npt, bmat, zmat, idz, ndim, vlag, beta,knew,w)
      call updateq(n, npt, idz, knew, prederr(1), xpt(knew, :),         &
     & bmat(knew, :), zmat, gq, hq, pq)

      ! Include the new interpolation point. This should be done after
      ! BMAT, ZMAT, and the model are updated.
      fval(knew) = f
      xpt(knew, :) = xnew

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
              itest = 0
          else
              gqsq = zero
              do i = 1, n
                  gqsq = gqsq + gq(i)**2
              end do
              vlag(1 : npt) = fval - fval(kopt)
               
!              galt = matmul(vlag(1 : npt), bmat(1 : npt, 1 : n))
!              galtsq = dot_product(galt, galt)
              galt = zero
              do k = 1, npt
                  galt = galt + vlag(k)*bmat(k, :)
              end do
              galtsq = zero
              do i = 1, n
                  galtsq = galtsq + galt(i)*galt(i) 
              end do

              ! Test whether to replace the new quadratic model by the
              ! least Frobenius norm interpolant, making the replacement
              ! if the test is satisfied.
              if (gqsq < 100.0_rp*galtsq) then
                  itest = 0
              else
                  itest = itest + 1
              end if

              if (itest >= 3) then
                  call qalt(gq, hq, pq, fval, bmat(1 : npt, :), zmat,   &
     &             n,npt,kopt,idz)
                  itest = 0
              end if
          end if
      end if

      if (f < fsave) kopt = knew
    
      ! If a trust region step has provided a sufficient decrease in F,
      ! then branch for another trust region calculation. The case
      ! KSAVE>0 occurs when the new function value was calculated by
      ! a model step.
      if (f <= fsave + tenth*vquad .or. ksav > 0) goto 100
      
      ! Alternatively, find out if the interpolation points are close
      ! enough to the best point so far.
      knew = 0
  460 distsq = 4.0_rp*delta*delta
      do k = 1, npt
          xdiff = xpt(k, :) - xopt
          xdsq(k) = dot_product(xdiff, xdiff)
      end do
      if (maxval(xdsq) > distsq) then
          knew = maxloc(xdsq, 1)
          distsq = maxval(xdsq)
      end if

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
          nfsav = nf
          
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
