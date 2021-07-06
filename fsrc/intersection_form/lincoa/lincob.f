!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of lincob.f90.
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


      subroutine lincob(n, npt, m, amat, b, x, rhobeg, rhoend, iprint, m&
     &axfun, xbase, xpt, fval, xsav, xopt, gopt, hq, pq, bmat, zmat, ndi&
     &m, step, sp, xnew, iact, rescon, qfac, rfac, pqw, w, f, info, ftar&
     &get)

! Dummy variables
      use consts_mod, only : IK, RP
      implicit none
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: m
      integer(IK), intent(in) :: n
      integer(IK), intent(in) :: ndim
      integer(IK), intent(in) :: npt
      integer(IK), intent(out) :: iact(:)
      integer(IK), intent(in) :: maxfun
      integer(IK), intent(out) :: info
      real(RP), intent(in) :: ftarget
      real(RP), intent(in) :: rhobeg
      real(RP), intent(out) :: rfac
      real(RP), intent(inout) :: amat(n, m)
      real(RP), intent(out) :: qfac(n, n)
      real(RP), dimension(npt, *) :: zmat
      real(RP), intent(in) :: rhoend
      real(RP), intent(inout) :: f
      real(RP), intent(inout), dimension(*) :: b
      real(RP), intent(inout), dimension(*) :: fval
      real(RP), intent(inout), dimension(*) :: gopt
      real(RP), intent(inout), dimension(*) :: hq
      real(RP), intent(inout), dimension(*) :: pq
      real(RP), intent(inout), dimension(*) :: pqw
      real(RP), intent(inout), dimension(*) :: rescon
      real(RP), intent(inout), dimension(*) :: sp
      real(RP), intent(inout), dimension(*) :: step
      real(RP), intent(inout), dimension(*) :: w
      real(RP), intent(inout), dimension(*) :: x
      real(RP), intent(inout), dimension(*) :: xbase
      real(RP), intent(inout), dimension(*) :: xnew
      real(RP), intent(inout), dimension(*) :: xopt
      real(RP), intent(inout), dimension(*) :: xsav
      real(RP), intent(inout), dimension(ndim, *) :: bmat
      real(RP), intent(inout), dimension(npt, *) :: xpt
! Local variables
      real(RP) :: almost_infinity, del, delsav, delta, dffalt, diff, dis&
     &tsq, fopt, fsave, half, one, qoptsq, ratio, rho, snorm, ssq, sum, &
     &sumz, temp, tenth, vqalt, vquad, xdiff, xoptsq, zero
      integer :: i, idz, ifeas, ih, imprv, ip, itest, j, k, knew, kopt, &
     &ksave, nact, nf, nh, np, nptm, nvala, nvalb
!*++
!*++ end of declarations rewritten by spag
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
!       identical to the corresponding arguments in SUBROUTINE LINCOA.
!     AMAT is a matrix whose columns are the constraint gradients, scaled
!       so that they have unit length.
!     B contains on entry the right hand sides of the constraints, scaled
!       as above, but later B is modified for variables relative to XBASE.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT contains the interpolation point coordinates relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XSAV holds the best feasible vector of variables so far, without any
!       shift of origin.
!     XOPT is set to XSAV-XBASE, which is the displacement from XBASE of
!       the feasible vector of variables that provides the least calculated
!       F so far, this vector being the current trust region centre.
!     GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of the big inverse matrix H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix
!       of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
!       where the elements of DZ are plus or minus one, as specified by IDZ.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     STEP is employed for trial steps from XOPT. It is also used for working
!       space when XBASE is shifted and in PRELIM.
!     SP is reserved for the scalar products XOPT^T XPT(K,.), K=1,2,...,NPT,
!       followed by STEP^T XPT(K,.), K=1,2,...,NPT.
!     XNEW is the displacement from XBASE of the vector of variables for
!       the current calculation of F, except that SUBROUTINE TRSTEP uses it
!       for working space.
!     IACT is an integer array for the indices of the active constraints.
!     RESCON holds useful information about the constraint residuals. Every
!       nonnegative RESCON(J) is the residual of the J-th constraint at the
!       current trust region centre. Otherwise, if RESCON(J) is negative, the
!       J-th constraint holds as a strict inequality at the trust region
!       centre, its residual being at least |RESCON(J)|; further, the value
!       of |RESCON(J)| is at least the current trust region radius DELTA.
!     QFAC is the orthogonal part of the QR factorization of the matrix of
!       active constraint gradients, these gradients being ordered in
!       accordance with IACT. When NACT is less than N, columns are added
!       to QFAC to complete an N by N orthogonal matrix, which is important
!       for keeping calculated steps sufficiently close to the boundaries
!       of the active constraints.
!     RFAC is the upper triangular part of this QR factorization, beginning
!       with the first diagonal element, followed by the two elements in the
!       upper triangular part of the second column and so on.
!     PQW is used for working space, mainly for storing second derivative
!       coefficients of quadratic functions. Its length is NPT+N.
!     The array W is also used for working space. The required number of
!       elements, namely MAX[M+3*N,2*M+N,2*NPT], is set in LINCOA.
!
!     Set some constants.
!
      half = 0.5D0
      one = 1.0D0
      tenth = 0.1D0
      zero = 0.0D0
      np = N + 1
      nh = (N * np) / 2
      nptm = Npt - np
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      almost_infinity = huge(0.0D0) / 2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 15-08-2019
! See the comments below line number 210
      imprv = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
!       ZMAT and SP for the first iteration. An important feature is that,
!       if the interpolation point XPT(K,.) is not feasible, where K is any
!       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
!       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
!       is set so that XPT(KOPT,.) is the initial trust region centre.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  STEP,PQW,W)
      call PRELIM(N, Npt, M, Amat, B, X, Rhobeg, Iprint, Xbase, Xpt, Fva&
     &l, Xsav, Xopt, Gopt, kopt, Hq, Pq, Bmat, Zmat, idz, Ndim, Sp, Resc&
     &on, Step, Pqw, W, F, Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 04-06-2019):
      if (F /= F .or. F > almost_infinity) then
          fopt = Fval(kopt)
          Info = -2
          goto 600
      end if
!     By Tom/Zaikun (on 04-06-2019/07-06-2019):
!     Note that we should NOT compare F and FTARGET, because X may not
!     be feasible at the exit of PRELIM.
      if (Fval(kopt) <= Ftarget) then
          F = Fval(kopt)
          X(1:N) = Xsav(1:N)
          Info = 1
          goto 700
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Begin the iterative procedure.
!
      nf = Npt
      fopt = Fval(kopt)
      rho = Rhobeg
      delta = rho
      ifeas = 0
      nact = 0
      itest = 3
100   knew = 0
      nvala = 0
      nvalb = 0
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!       to BMAT that do not depend on ZMAT.
!
200   fsave = fopt
      xoptsq = zero
      do i = 1, N
          xoptsq = xoptsq + Xopt(i)**2
      end do
      if (xoptsq >= 1.0D4 * delta * delta) then
          qoptsq = 0.25D0 * xoptsq
          do k = 1, Npt
              sum = zero
              do i = 1, N
                  sum = sum + Xpt(k, i) * Xopt(i)
              end do
              sum = sum - half * xoptsq
              W(Npt + k) = sum
              Sp(k) = zero
              do i = 1, N
                  Xpt(k, i) = Xpt(k, i) - half * Xopt(i)
                  Step(i) = Bmat(k, i)
                  W(i) = sum * Xpt(k, i) + qoptsq * Xopt(i)
                  ip = Npt + i
                  do j = 1, i
                      Bmat(ip, j) = Bmat(ip, j) + Step(i) * W(j) + W(i) &
     &* Step(j)
                  end do
              end do
          end do
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
          do k = 1, nptm
              sumz = zero
              do i = 1, Npt
                  sumz = sumz + Zmat(i, k)
                  W(i) = W(Npt + i) * Zmat(i, k)
              end do
              do j = 1, N
                  sum = qoptsq * sumz * Xopt(j)
                  do i = 1, Npt
                      sum = sum + W(i) * Xpt(i, j)
                  end do
                  Step(j) = sum
                  if (k < idz) sum = -sum
                  do i = 1, Npt
                      Bmat(i, j) = Bmat(i, j) + sum * Zmat(i, k)
                  end do
              end do
              do i = 1, N
                  ip = i + Npt
                  temp = Step(i)
                  if (k < idz) temp = -temp
                  do j = 1, i
                      Bmat(ip, j) = Bmat(ip, j) + temp * Step(j)
                  end do
              end do
          end do
!
!     Update the right hand sides of the constraints.
!
          if (M > 0) then
              do j = 1, M
                  temp = zero
                  do i = 1, N
                      temp = temp + Amat(i, j) * Xopt(i)
                  end do
                  B(j) = B(j) - temp
              end do
          end if
!
!     The following instructions complete the shift of XBASE, including the
!       changes to the parameters of the quadratic model.
!
          ih = 0
          do j = 1, N
              W(j) = zero
              do k = 1, Npt
                  W(j) = W(j) + Pq(k) * Xpt(k, j)
                  Xpt(k, j) = Xpt(k, j) - half * Xopt(j)
              end do
              do i = 1, j
                  ih = ih + 1
                  Hq(ih) = Hq(ih) + W(i) * Xopt(j) + Xopt(i) * W(j)
                  Bmat(Npt + i, j) = Bmat(Npt + j, i)
              end do
          end do
          do j = 1, N
              Xbase(j) = Xbase(j) + Xopt(j)
              Xopt(j) = zero
              Xpt(kopt, j) = zero
          end do
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 21-03-2020
! Exit if BMAT or ZMAT contians NaN
      do j = 1, N
          do i = 1, Ndim
              if (Bmat(i, j) /= Bmat(i, j)) then
                  Info = -3
                  goto 600
              end if
          end do
      end do
      do j = 1, nptm
          do i = 1, Npt
              if (Zmat(i, j) /= Zmat(i, j)) then
                  Info = -3
                  goto 600
              end if
          end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     In the case KNEW=0, generate the next trust region step by calling
!       TRSTEP, where SNORM is the current trust region radius initially.
!       The final value of SNORM is the length of the calculated step,
!       except that SNORM is zero on return if the projected gradient is
!       unsuitable for starting the conjugate gradient iterations.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRSTEP or QMSTEP is not predictable, and Segmentation Fault or
! infinite cycling may happen. This is because any equality/inequality
! comparison involving NaN returns FALSE, which can lead to unintended
! behavior of the code, including uninitialized indices, which can lead
! to segmentation faults.
      do j = 1, N
          if (Gopt(j) /= Gopt(j)) then
              Info = -3
              goto 600
          end if
      end do
      do i = 1, nh
          if (Hq(i) /= Hq(i)) then
              Info = -3
              goto 600
          end if
      end do
      do i = 1, Npt
          if (Pq(i) /= Pq(i)) then
              Info = -3
              goto 600
          end if
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      delsav = delta
      ksave = knew
      if (knew == 0) then
          snorm = delta
          do i = 1, N
              Xnew(i) = Gopt(i)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 19-03-2020: B is never used in TRSTEP
!          CALL TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call TRSTEP(N, Npt, M, Amat, Xpt, Hq, Pq, nact, Iact, Rescon, &
     &Qfac, Rfac, snorm, Step, Xnew, W, W(M + 1), Pqw, Pqw(np), W(M + np&
     &))
!
!     A trust region step is applied whenever its length, namely SNORM, is at
!       least HALF*DELTA. It is also applied if its length is at least 0.1999
!       times DELTA and if a line search of TRSTEP has caused a change to the
!       active set. Otherwise there is a branch below to label 530 or 560.
!
          temp = half * delta
          if (Xnew(1) >= half) temp = 0.1999D0 * delta
          if (snorm <= temp) then
              delta = half * delta
              if (delta <= 1.4D0 * rho) delta = rho
              nvala = nvala + 1
              nvalb = nvalb + 1
              temp = snorm / rho
              if (delsav > rho) temp = one
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 24-07-2019
!              IF (TEMP .GE. HALF) NVALA=ZERO
!              IF (TEMP .GE. TENTH) NVALB=ZERO
              if (temp >= half) nvala = 0
              if (temp >= tenth) nvalb = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (delsav > rho) goto 400
              if (nvala < 5 .and. nvalb < 3) goto 400
              if (snorm > zero) ksave = -1
              goto 500
          end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 24-07-2019
!          NVALA=ZERO
!          NVALB=ZERO
          nvala = 0
          nvalb = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Alternatively, KNEW is positive. Then the model step is calculated
!       within a trust region of radius DEL, after setting the gradient at
!       XBASE and the second derivative parameters of the KNEW-th Lagrange
!       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
!
      else
          del = DMAX1(tenth * delta, rho)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: See the comments below line number 140
!          DO 160 I=1,N
!  160     W(I)=BMAT(KNEW,I)
          do i = 1, N
              W(i) = Bmat(knew, i)
              if (W(i) /= W(i)) then
                  Info = -3
                  goto 600
              end if
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k = 1, Npt
              Pqw(k) = zero
          end do
          do j = 1, nptm
              temp = Zmat(knew, j)
              if (j < idz) temp = -temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: See the comments below line number 140
! Note that the data in PQW is used in QMSTEP below
!          DO 180 K=1,NPT
!  180     PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
              do k = 1, Npt
                  Pqw(k) = Pqw(k) + temp * Zmat(k, j)
                  if (Pqw(k) /= Pqw(k)) then
                      Info = -3
                      goto 600
                  end if
              end do
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: B is never used in QMSTEP
!          CALL QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call QMSTEP(N, Npt, M, Amat, Xpt, Xopt, nact, Iact, Rescon, Qf&
     &ac, kopt, knew, del, Step, W, Pqw, W(np), W(np + M), ifeas)
      end if
!
!     Set VQUAD to the change to the quadratic model when the move STEP is
!       made from XOPT. If STEP is a trust region step, then VQUAD should be
!       negative. If it is nonnegative due to rounding errors in this case,
!       there is a branch to label 530 to try to improve the model.
!
      vquad = zero
      ih = 0
      do j = 1, N
          vquad = vquad + Step(j) * Gopt(j)
          do i = 1, j
              ih = ih + 1
              temp = Step(i) * Step(j)
              if (i == j) temp = half * temp
              vquad = vquad + temp * Hq(ih)
          end do
      end do
      do k = 1, Npt
          temp = zero
          do j = 1, N
              temp = temp + Xpt(k, j) * Step(j)
              Sp(Npt + k) = temp
          end do
          vquad = vquad + half * Pq(k) * temp * temp
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 15-08-2019
! Although very rarely, with the original code, an infinite loop can occur
! in the following scenario.
! Suppose that, at an certain iteration,
! KNEW = 0, SNORM > 0.5*DELTA > RHO, VQUAD >= 0, and
! sum_{K=1}^NPT ||XPT(K,:)-XOPT(:)||^2 < DELTA^2
! (i.e., DELTA is large and SNORM is not small, yet VQUAD >= 0 due to
! rounding errors and XPT are not far from XOPT).
! Then the program will goto 530 and then goto 20, where XBASE may be
! shifted to the current best point, in the hope of reducing rounding
! errors and 'improve' the model. Afterwards, another trust region step
! is produced by the 'improved' model. Note that DELTA remains unchanged
! in this process. If the new trust region step turns out to satisfy
! SNORM > 0.5*DELTA and VQUAD >= 0 again (i.e., the 'improved' model
! still suffers from rounding errors), then the program will goto 530
! and then goto 20, where shifting will not happen because either XBASE
! was already shifted to the current best point in last step, or XBASE
! is close to the current best point. Consequently, the model will
! remain unchanged, and produce the same trust region step again. This
! leads to an infinite loop.
! The infinite loop did happen when the MATLAB interface was applied to
! min atan(x+100) s.t. x<=-99 (x0=-99, npt=3, rhobeg=1, rhoend=1e-6).
! The problem does not exist in NEWUOA or BOBYQA, where the program will
! exit immediately when VQUAD >= 0.
! To prevent such a loop, here we use IMPRV to record whether the path
! 530 --> 20 has already happened for last trust region step. IMPRV=1
! implies that last trust region step satisfies VQUAD >= 0 and followed
! 530 --> 20. With IMPRV=1, if VQUAD is again nonnegative for the new trust
! region step, we should not goto 530 but goto 560, where IMPRV will be
! set to 0 and DELTA will be reduced. Otherwise, an infinite loop would happen.
!      IF (KSAVE .EQ. 0 .AND. VQUAD .GE. ZERO) GOTO 530
      if (ksave == 0 .and. .not. (vquad < zero)) then
          if (imprv == 1) goto 500
          imprv = 1
          goto 400
      else
          imprv = 0
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the next value of the objective function. The difference
!       between the actual new value of F and the value predicted by the
!       model is recorded in DIFF.
!
300   nf = nf + 1
      if (nf > Maxfun) then
          nf = nf - 1
          if (Iprint > 0) print 99001
99001 format(/4X, 'Return from LINCOA because CALFUN has been', ' called&
     & MAXFUN times.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          Info = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          goto 600
      end if
      xdiff = zero
      do i = 1, N
          Xnew(i) = Xopt(i) + Step(i)
          X(i) = Xbase(i) + Xnew(i)
          xdiff = xdiff + (X(i) - Xsav(i))**2
      end do
      xdiff = DSQRT(xdiff)
      if (ksave == -1) xdiff = rho
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (XDIFF .LE. TENTH*RHO .OR. XDIFF .GE. DELTA+DELTA) THEN
      if (.not. (xdiff > tenth * rho .and. xdiff < delta + delta)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ifeas = 0
          if (Iprint > 0) print 99002
99002 format(/4X, 'Return from LINCOA because rounding errors', ' preven&
     &t reasonable changes to X.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          Info = 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          goto 600
      end if
      if (ksave <= 0) ifeas = 1
      F = DFLOAT(ifeas)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i = 1, N
          if (X(i) /= X(i)) then
              F = X(i)
              ! Set F to NaN
              if (nf == 1) then
                  fopt = F
                  Xopt(1:N) = zero
              end if
              Info = -1
              goto 600
          end if
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call CALFUN(N, X, F)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 04-06-2019):
      if (F /= F .or. F > almost_infinity) then
          if (nf == 1) then
              fopt = F
              Xopt(1:N) = zero
          end if
          Info = -2
          goto 600
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (Iprint == 3) then
          print 99003, nf, F, (X(i), i=1, N)
99003 format(/4X, 'Function number', I6, ' F =', 1PD18.10, ' The corresp&
     &onding X is:'/(2X, 5D15.6))
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (KSAVE .EQ. -1) GOTO 600
      if (ksave == -1) then
          Info = 0
          goto 600
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      diff = F - fopt - vquad
!
!     If X is feasible, then set DFFALT to the difference between the new
!       value of F and the value predicted by the alternative model.
!
      if (ifeas == 1 .and. itest < 3) then
          do k = 1, Npt
              Pqw(k) = zero
              W(k) = Fval(k) - Fval(kopt)
          end do
          do j = 1, nptm
              sum = zero
              do i = 1, Npt
                  sum = sum + W(i) * Zmat(i, j)
              end do
              if (j < idz) sum = -sum
              do k = 1, Npt
                  Pqw(k) = Pqw(k) + sum * Zmat(k, j)
              end do
          end do
          vqalt = zero
          do k = 1, Npt
              sum = zero
              do j = 1, N
                  sum = sum + Bmat(k, j) * Step(j)
              end do
              vqalt = vqalt + sum * W(k)
              vqalt = vqalt + Pqw(k) * Sp(Npt + k) * (half * Sp(Npt + k)&
     & + Sp(k))
          end do
          dffalt = F - fopt - vqalt
      end if
      if (itest == 3) then
          dffalt = diff
          itest = 0
      end if
!
!     Pick the next value of DELTA after a trust region step.
!
      if (ksave == 0) then
          ratio = (F - fopt) / vquad
          if (ratio <= tenth) then
              delta = half * delta
          elseif (ratio <= 0.7D0) then
              delta = DMAX1(half * delta, snorm)
          else
              temp = DSQRT(2.0D0) * delta
              delta = DMAX1(half * delta, snorm + snorm)
              delta = DMIN1(delta, temp)
          end if
          if (delta <= 1.4D0 * rho) delta = rho
      end if
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!       can be moved. If STEP is a trust region step, then KNEW is zero at
!       present, but a positive value is picked by subroutine UPDATE.
!
      call UPDATE(N, Npt, Xpt, Bmat, Zmat, idz, Ndim, Sp, Step, kopt, kn&
     &ew, Pqw, W)
      if (knew == 0) then
          if (Iprint > 0) print 99004
99004 format(/4X, 'Return from LINCOA because the denominator'' of the u&
     &pdating form ula is zero.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          Info = 9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          goto 600
      end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 19-03-2020
! Exit if BMAT or ZMAT contians NaN
      do j = 1, N
          do i = 1, Ndim
              if (Bmat(i, j) /= Bmat(i, j)) then
                  Info = -3
                  goto 600
              end if
          end do
      end do
      do j = 1, nptm
          do i = 1, Npt
              if (Zmat(i, j) /= Zmat(i, j)) then
                  Info = -3
                  goto 600
              end if
          end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     If ITEST is increased to 3, then the next quadratic model is the
!       one whose second derivative matrix is least subject to the new
!       interpolation conditions. Otherwise the new model is constructed
!       by the symmetric Broyden method in the usual way.
!
      if (ifeas == 1) then
          itest = itest + 1
          if (DABS(dffalt) >= tenth * DABS(diff)) itest = 0
      end if
!
!     Update the second derivatives of the model by the symmetric Broyden
!       method, using PQW for the second derivative parameters of the new
!       KNEW-th Lagrange function. The contribution from the old parameter
!       PQ(KNEW) is included in the second derivative matrix HQ. W is used
!       later for the gradient of the new KNEW-th Lagrange function.
!
      if (itest < 3) then
          do k = 1, Npt
              Pqw(k) = zero
          end do
          do j = 1, nptm
              temp = Zmat(knew, j)
              if (temp /= zero) then
                  if (j < idz) temp = -temp
                  do k = 1, Npt
                      Pqw(k) = Pqw(k) + temp * Zmat(k, j)
                  end do
              end if
          end do
          ih = 0
          do i = 1, N
              W(i) = Bmat(knew, i)
              temp = Pq(knew) * Xpt(knew, i)
              do j = 1, i
                  ih = ih + 1
                  Hq(ih) = Hq(ih) + temp * Xpt(knew, j)
              end do
          end do
          Pq(knew) = zero
          do k = 1, Npt
              Pq(k) = Pq(k) + diff * Pqw(k)
          end do
      end if
!
!     Include the new interpolation point with the corresponding updates of
!       SP. Also make the changes of the symmetric Broyden method to GOPT at
!       the old XOPT if ITEST is less than 3.
!
      Fval(knew) = F
      Sp(knew) = Sp(kopt) + Sp(Npt + kopt)
      ssq = zero
      do i = 1, N
          Xpt(knew, i) = Xnew(i)
          ssq = ssq + Step(i)**2
      end do
      Sp(Npt + knew) = Sp(Npt + kopt) + ssq
      if (itest < 3) then
          do k = 1, Npt
              temp = Pqw(k) * Sp(k)
              do i = 1, N
                  W(i) = W(i) + temp * Xpt(k, i)
              end do
          end do
          do i = 1, N
              Gopt(i) = Gopt(i) + diff * W(i)
          end do
      end if
!
!     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
!       least calculated value so far with a feasible vector of variables.
!
      if (F < fopt .and. ifeas == 1) then
          fopt = F
          do j = 1, N
              Xsav(j) = X(j)
              Xopt(j) = Xnew(j)
          end do
          kopt = knew
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!         By Tom (on 04-06-2019):
          if (fopt <= Ftarget) then
              Info = 1
              goto 700
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          snorm = DSQRT(ssq)
          do j = 1, M
              if (Rescon(j) >= delta + snorm) then
                  Rescon(j) = snorm - Rescon(j)
              else
                  Rescon(j) = Rescon(j) + snorm
                  if (Rescon(j) + delta > zero) then
                      temp = B(j)
                      do i = 1, N
                          temp = temp - Xopt(i) * Amat(i, j)
                      end do
                      temp = DMAX1(temp, zero)
                      if (temp >= delta) temp = -temp
                      Rescon(j) = temp
                  end if
              end if
          end do
          do k = 1, Npt
              Sp(k) = Sp(k) + Sp(Npt + k)
          end do
!
!     Also revise GOPT when symmetric Broyden updating is applied.
!
          if (itest < 3) then
              ih = 0
              do j = 1, N
                  do i = 1, j
                      ih = ih + 1
                      if (i < j) Gopt(j) = Gopt(j) + Hq(ih) * Step(i)
                      Gopt(i) = Gopt(i) + Hq(ih) * Step(j)
                  end do
              end do
              do k = 1, Npt
                  temp = Pq(k) * Sp(Npt + k)
                  do i = 1, N
                      Gopt(i) = Gopt(i) + temp * Xpt(k, i)
                  end do
              end do
          end if
      end if
!
!     Replace the current model by the least Frobenius norm interpolant if
!       this interpolant gives substantial reductions in the predictions
!       of values of F at feasible points.
!
      if (itest == 3) then
          do k = 1, Npt
              Pq(k) = zero
              W(k) = Fval(k) - Fval(kopt)
          end do
          do j = 1, nptm
              sum = zero
              do i = 1, Npt
                  sum = sum + W(i) * Zmat(i, j)
              end do
              if (j < idz) sum = -sum
              do k = 1, Npt
                  Pq(k) = Pq(k) + sum * Zmat(k, j)
              end do
          end do
          do j = 1, N
              Gopt(j) = zero
              do i = 1, Npt
                  Gopt(j) = Gopt(j) + W(i) * Bmat(i, j)
              end do
          end do
          do k = 1, Npt
              temp = Pq(k) * Sp(k)
              do i = 1, N
                  Gopt(i) = Gopt(i) + temp * Xpt(k, i)
              end do
          end do
          do ih = 1, nh
              Hq(ih) = zero
          end do
      end if
!
!     If a trust region step has provided a sufficient decrease in F, then
!       branch for another trust region calculation. Every iteration that
!       takes a model step is followed by an attempt to take a trust region
!       step.
!
      knew = 0
      if (ksave > 0) goto 200
      if (ratio >= tenth) goto 200
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
400   distsq = DMAX1(delta * delta, 4.0D0 * rho * rho)
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
!     If KNEW is positive, then branch back for the next iteration, which
!       will generate a "model step". Otherwise, if the current iteration
!       has reduced F, or if DELTA was above its lower bound when the last
!       trust region step was calculated, then try a "trust region" step
!       instead.
!
      if (knew > 0) goto 200
      knew = 0
      if (fopt < fsave) goto 200
      if (delsav > rho) goto 200
!
!     The calculations with the current value of RHO are complete.
!       Pick the next value of RHO.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 15-08-2019
! See the comments below line number 210
!  560 IF (RHO .GT. RHOEND) THEN
500   imprv = 0
      if (rho > Rhoend) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          delta = half * rho
          if (rho > 250.0D0 * Rhoend) then
              rho = tenth * rho
          elseif (rho <= 16.0D0 * Rhoend) then
              rho = Rhoend
          else
              rho = DSQRT(rho * Rhoend)
          end if
          delta = DMAX1(delta, rho)
          if (Iprint >= 2) then
              if (Iprint >= 3) print 99005
99005 format(5X)
              print 99006, rho, nf
99006 format(/4X, 'New RHO =', 1PD11.4, 5X, 'Number of', ' function valu&
     &es =', I6)
              print 99008, fopt, (Xbase(i) + Xopt(i), i=1, N)
          end if
          goto 100
      end if
!
!     Return from the calculation, after branching to label 220 for another
!       Newton-Raphson step if it has not been tried before.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ksave == -1) goto 300
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  600 IF (FOPT .LE. F .OR. IFEAS .EQ. 0) THEN
600   if (fopt <= F .or. ifeas == 0 .or. F /= F) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i = 1, N
              X(i) = Xsav(i)
          end do
          F = fopt
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (IPRINT .GE. 1) THEN
700   if (Iprint >= 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          print 99007, nf
99007 format(/4X, 'At the return from LINCOA', 5X, 'Number of function v&
     &alues =', I6)
          print 99008, F, (X(i), i=1, N)
      end if
      W(1) = F
      W(2) = DFLOAT(nf) + half
99008 format(4X, 'Least value of F =', 1PD23.15, 9X, 'The corresponding &
     &X is:'/(2X, 5D15.6))
      end subroutine LINCOB