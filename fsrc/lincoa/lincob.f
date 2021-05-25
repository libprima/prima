      subroutine LINCOB(N, NPT, M, AMAT, B, X, RHOBEG, RHOEND, IPRINT,
1     MAXFUN, XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM,
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C 2 STEP, SP, XNEW, IACT, RESCON, QFAC, RFAC, PQW, W)
2     STEP, SP, XNEW, IACT, RESCON, QFAC, RFAC, PQW, W, F, INFO, FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C implicit real * 8(A - H, O - Z)
      implicit real(kind(0.0D0)) (A - H, O - Z)
      implicit integer(I - N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension AMAT(N, *), B(*), X(*), XBASE(*), XPT(NPT, *), FVAL(*),
1     XSAV(*), XOPT(*), GOPT(*), HQ(*), PQ(*), BMAT(NDIM, *),
2     ZMAT(NPT, *), STEP(*), SP(*), XNEW(*), IACT(*), RESCON(*),
3     QFAC(N, *), RFAC(*), PQW(*), W(*)
      C
      C The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
      C identical to the corresponding arguments in subroutine LINCOA.
      C AMAT is a matrix whose columns are the constraint gradients, scaled
      C so that they have unit length.
      C B contains on entry the right hand sides of the constraints, scaled
      C as above, but later B is modified for variables relative to XBASE.
      C XBASE holds a shift of origin that should reduce the contributions
      C from rounding errors to values of the model and Lagrange functions.
      C XPT contains the interpolation point coordinates relative to XBASE.
      C FVAL holds the values of F at the interpolation points.
      C XSAV holds the best feasible vector of variables so far, without any
      C shift of origin.
      C XOPT is set to XSAV - XBASE, which is the displacement from XBASE of
      C the feasible vector of variables that provides the least calculated
      C F so far, this vector being the current trust region centre.
      C GOPT holds the gradient of the quadratic model at XSAV = XBASE + XOPT.
      C HQ holds the explicit second derivatives of the quadratic model.
      C PQ contains the parameters of the implicit second derivatives of the
      C quadratic model.
      C BMAT holds the last N columns of the big inverse matrix H.
      C ZMAT holds the factorization of the leading NPT by NPT submatrix
      C of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
      C where the elements of DZ are plus or minus one, as specified by IDZ.
      C NDIM is the first dimension of BMAT and has the value NPT + N.
      C STEP is employed for trial steps from XOPT.It is also used for working
      C space when XBASE is shifted and in PRELIM.
      C SP is reserved for the scalar products XOPT^T XPT(K, .), K = 1, 2, ..., NPT,
      C followed by STEP^T XPT(K, .), K = 1, 2, ..., NPT.
      C XNEW is the displacement from XBASE of the vector of variables for
      C the current calculation of F, except that subroutine TRSTEP uses it
      C for working space.
      C IACT is an integer array for the indices of the active constraints.
      C RESCON holds useful information about the constraint residuals.Every
      C nonnegative RESCON(J) is the residual of the J - th constraint at the
      C current trust region centre.Otherwise, if RESCON(J) is negative, the
      C J - th constraint holds as a strict inequality at the trust region
      C centre, its residual being at least|RESCON(J) |; further, the value
      C of|RESCON(J) |is at least the current trust region radius DELTA.
      C QFAC is the orthogonal part of the QR factorization of the matrix of
      C active constraint gradients, these gradients being ordered in
      C accordance with IACT.When NACT is less than N, columns are added
      C to QFAC to complete an N by N orthogonal matrix, which is important
      C for keeping calculated steps sufficiently close to the boundaries
      C of the active constraints.
      C RFAC is the upper triangular part of this QR factorization, beginning
      C with the first diagonal element, followed by the two elements in the
      C upper triangular part of the second column and so on.
      C PQW is used for working space, mainly for storing second derivative
      C coefficients of quadratic functions.Its length is NPT + N.
      C The array W is also used for working space.The required number of
      C elements, namely MAX[M + 3 * N, 2 * M + N, 2 * NPT], is set in LINCOA.
      C
      C Set some constants.
      C
      HALF = 0.5D0
      ONE = 1.0D0
      TENTH = 0.1D0
      ZERO = 0.0D0
      NP = N + 1
      NH = (N * NP) / 2
      NPTM = NPT - NP
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      ALMOST_INFINITY = huge(0.0D0) / 2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C Zaikun 15 - 08 - 2019
      C See the comments below line number 210
      IMPRV = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      C
      C Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
      C ZMAT and SP for the first iteration.An important feature is that,
      C if the interpolation point XPT(K, .) is not feasible, where K is any
      C integer from[1, NPT], then a change is made to XPT(K, .) if necessary
      C so that the constraint violation is at least 0.2 * RHOBEG.Also KOPT
      C is set so that XPT(KOPT, .) is the initial trust region centre.
      C
      call PRELIM(N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
1     XSAV, XOPT, GOPT, KOPT, HQ, PQ, BMAT, ZMAT, IDZ, NDIM, SP, RESCON,
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C 2 STEP, PQW, W)
2     STEP, PQW, W, F, FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C By Tom(on 04 - 06 - 2019):
      if (F /= F .or. F > ALMOST_INFINITY) then
          FOPT = FVAL(KOPT)
          INFO = -2
          goto 600
      end if
      C By Tom / Zaikun(on 04 - 06 - 2019 / 07 - 06 - 2019):
      C Note that we should NOT compare F and FTARGET, because X may not
      C be feasible at the exit of PRELIM.
      if (FVAL(KOPT) <= FTARGET) then
          F = FVAL(KOPT)
          X(1:N) = XSAV(1:N)
          INFO = 1
          goto 616
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      C
      C Begin the iterative procedure.
      C
      NF = NPT
      FOPT = FVAL(KOPT)
      RHO = RHOBEG
      DELTA = RHO
      IFEAS = 0
      NACT = 0
      ITEST = 3
10    KNEW = 0
      NVALA = 0
      NVALB = 0
      C
      C Shift XBASE if XOPT may be too far from XBASE.First make the changes
      C to BMAT that do not depend on ZMAT.
      C
20    FSAVE = FOPT
      XOPTSQ = ZERO
      do I = 1, N
          XOPTSQ = XOPTSQ + XOPT(I)**2
      end do
      if (XOPTSQ >= 1.0D4 * DELTA * DELTA) then
          QOPTSQ = 0.25D0 * XOPTSQ
          do K = 1, NPT
              SUM = ZERO
              do I = 1, N
                  SUM = SUM + XPT(K, I) * XOPT(I)
              end do
              SUM = SUM - HALF * XOPTSQ
              W(NPT + K) = SUM
              SP(K) = ZERO
              do I = 1, N
                  XPT(K, I) = XPT(K, I) - HALF * XOPT(I)
                  STEP(I) = BMAT(K, I)
                  W(I) = SUM * XPT(K, I) + QOPTSQ * XOPT(I)
                  IP = NPT + I
                  do J = 1, I
                      BMAT(IP, J) = BMAT(IP, J) + STEP(I) * W(J) + W(I) * STEP(J)
                  end do
              end do
          end do
          C
          C then the revisions of BMAT that depend on ZMAT are calculated.
          C
          do K = 1, NPTM
              SUMZ = ZERO
              do I = 1, NPT
                  SUMZ = SUMZ + ZMAT(I, K)
                  W(I) = W(NPT + I) * ZMAT(I, K)
              end do
              do J = 1, N
                  SUM = QOPTSQ * SUMZ * XOPT(J)
                  do I = 1, NPT
                      SUM = SUM + W(I) * XPT(I, J)
                  end do
                  STEP(J) = SUM
                  if (K < IDZ) SUM = -SUM
                  do I = 1, NPT
                      BMAT(I, J) = BMAT(I, J) + SUM * ZMAT(I, K)
                  end do
              end do
              do I = 1, N
                  IP = I + NPT
                  TEMP = STEP(I)
                  if (K < IDZ) TEMP = -TEMP
                  do J = 1, I
                      BMAT(IP, J) = BMAT(IP, J) + TEMP * STEP(J)
                  end do
              end do
          end do
          C
          C Update the right hand sides of the constraints.
          C
          if (M > 0) then
              do J = 1, M
                  TEMP = ZERO
                  do I = 1, N
                      TEMP = TEMP + AMAT(I, J) * XOPT(I)
                  end do
                  B(J) = B(J) - TEMP
              end do
          end if
          C
          C The following instructions complete the shift of XBASE, including the
          C changes to the parameters of the quadratic model.
          C
          IH = 0
          do J = 1, N
              W(J) = ZERO
              do K = 1, NPT
                  W(J) = W(J) + PQ(K) * XPT(K, J)
                  XPT(K, J) = XPT(K, J) - HALF * XOPT(J)
              end do
              do I = 1, J
                  IH = IH + 1
                  HQ(IH) = HQ(IH) + W(I) * XOPT(J) + XOPT(I) * W(J)
                  BMAT(NPT + I, J) = BMAT(NPT + J, I)
              end do
          end do
          do J = 1, N
              XBASE(J) = XBASE(J) + XOPT(J)
              XOPT(J) = ZERO
              XPT(KOPT, J) = ZERO
          end do
      end if

      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C Zaikun 21 - 03 - 2020
      C exit if BMAT or ZMAT contians NaN
      do J = 1, N
          do I = 1, NDIM
              if (BMAT(I, J) /= BMAT(I, J)) then
                  INFO = -3
                  goto 600
              end if
          end do
      end do
      do J = 1, NPTM
          do I = 1, NPT
              if (ZMAT(I, J) /= ZMAT(I, J)) then
                  INFO = -3
                  goto 600
              end if
          end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      C
      C In the case KNEW = 0, generate the next trust region step by calling
      C TRSTEP, where SNORM is the current trust region radius initially.
      C The final value of SNORM is the length of the calculated step,
      C except that SNORM is zero on return if the projected gradient is
      C unsuitable for starting the conjugate gradient iterations.
      C
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C Zaikun 2019 - 08 - 29:For ill - conditioned problems, NaN may occur in the
      C models.In such a case, we terminate the code.Otherwise, the behavior
      C of TRSTEP or QMSTEP is not predictable, and Segmentation Fault or
      C infinite cycling may happen.This is because any equality / inequality
      C comparison involving NaN returns FALSE, which can lead to unintended
      C behavior of the code, including uninitialized indices, which can lead
      C to segmentation faults.
      do J = 1, N
          if (GOPT(J) /= GOPT(J)) then
              INFO = -3
              goto 600
          end if
      end do
      do I = 1, NH
          if (HQ(I) /= HQ(I)) then
              INFO = -3
              goto 600
          end if
      end do
      do I = 1, NPT
          if (PQ(I) /= PQ(I)) then
              INFO = -3
              goto 600
          end if
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DELSAV = DELTA
      KSAVE = KNEW
      if (KNEW == 0) then
          SNORM = DELTA
          do I = 1, N
              XNEW(I) = GOPT(I)
          end do
          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          C Zaikun 19 - 03 - 2020:B is never used in TRSTEP
          C call TRSTEP(N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON,
          call TRSTEP(N, NPT, M, AMAT, XPT, HQ, PQ, NACT, IACT, RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1         QFAC, RFAC, SNORM, STEP, XNEW, W, W(M + 1), PQW, PQW(NP), W(M + NP))
          C
          C A trust region step is applied whenever its length, namely SNORM, is at
          C least HALF * DELTA.It is also applied if its length is at least 0.1999
          C times DELTA and if a line search of TRSTEP has caused a change to the
          C active set.Otherwise there is a branch below to label 530 or 560.
          C
          TEMP = HALF * DELTA
          if (XNEW(1) >= HALF) TEMP = 0.1999D0 * DELTA
          if (SNORM <= TEMP) then
              DELTA = HALF * DELTA
              if (DELTA <= 1.4D0 * RHO) DELTA = RHO
              NVALA = NVALA + 1
              NVALB = NVALB + 1
              TEMP = SNORM / RHO
              if (DELSAV > RHO) TEMP = ONE
              CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              C Zaikun 24 - 07 - 2019
              C if(TEMP >= HALF) NVALA = ZERO
              C if(TEMP >= TENTH) NVALB = ZERO
              if (TEMP >= HALF) NVALA = 0
              if (TEMP >= TENTH) NVALB = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (DELSAV > RHO) goto 530
              if (NVALA < 5 .and. NVALB < 3) goto 530
              if (SNORM > ZERO) KSAVE = -1
              goto 560
          end if
          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          C Zaikun 24 - 07 - 2019
          C NVALA = ZERO
          C NVALB = ZERO
          NVALA = 0
          NVALB = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          C
          C Alternatively, KNEW is positive.then the model step is calculated
          C within a trust region of radius DEL, after setting the gradient at
          C XBASE and the second derivative parameters of the KNEW - th Lagrange
          C function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
          C
      else
          DEL = DMAX1(TENTH * DELTA, RHO)
          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          C Zaikun 2019 - 08 - 29:See the comments below line number 140
          C do 160 I = 1, N
          C 160 W(I) = BMAT(KNEW, I)
          do I = 1, N
              W(I) = BMAT(KNEW, I)
              if (W(I) /= W(I)) then
                  INFO = -3
                  goto 600
              end if
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do K = 1, NPT
              PQW(K) = ZERO
          end do
          do J = 1, NPTM
              TEMP = ZMAT(KNEW, J)
              if (J < IDZ) TEMP = -TEMP
              CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              C Zaikun 2019 - 08 - 29:See the comments below line number 140
              C Note that the data in PQW is used in QMSTEP below
              C do 180 K = 1, NPT
              C 180 PQW(K) = PQW(K) + TEMP * ZMAT(K, J)
              do K = 1, NPT
                  PQW(K) = PQW(K) + TEMP * ZMAT(K, J)
                  if (PQW(K) /= PQW(K)) then
                      INFO = -3
                      goto 600
                  end if
              end do
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          C Zaikun 2019 - 08 - 29:B is never used in QMSTEP
          C call QMSTEP(N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON,
          call QMSTEP(N, NPT, M, AMAT, XPT, XOPT, NACT, IACT, RESCON,
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1         QFAC, KOPT, KNEW, DEL, STEP, W, PQW, W(NP), W(NP + M), IFEAS)
      end if
      C
      C Set VQUAD to the change to the quadratic model when the move STEP is
      C made from XOPT.if STEP is a trust region step, then VQUAD should be
      C negative.if it is nonnegative due to rounding errors in this case,
      C there is a branch to label 530 to try to improve the model.
      C
      VQUAD = ZERO
      IH = 0
      do J = 1, N
          VQUAD = VQUAD + STEP(J) * GOPT(J)
          do I = 1, J
              IH = IH + 1
              TEMP = STEP(I) * STEP(J)
              if (I == J) TEMP = HALF * TEMP
              VQUAD = VQUAD + TEMP * HQ(IH)
          end do
      end do
      do K = 1, NPT
          TEMP = ZERO
          do J = 1, N
              TEMP = TEMP + XPT(K, J) * STEP(J)
              SP(NPT + K) = TEMP
          end do
          VQUAD = VQUAD + HALF * PQ(K) * TEMP * TEMP
      end do
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C Zaikun 15 - 08 - 2019
      C Although very rarely, with the original code, an infinite loop can occur
      C in the following scenario.
      C Suppose that, at an certain iteration,
      C KNEW = 0, SNORM > 0.5 * DELTA > RHO, VQUAD >= 0, and
      C sum_{K = 1}^NPT||XPT(K, :) - XOPT(:) ||^2 < DELTA^2
      C(i.e., DELTA is large and SNORM is not small, yet VQUAD >= 0 due to
      C rounding errors and XPT are not far from XOPT) .
      C then the program will goto 530 and then goto 20, where XBASE may be
      C shifted to the current best point, in the hope of reducing rounding
      C errors and 'improve'the model.Afterwards, another trust region step
      C is produced by the 'improved'model.Note that DELTA remains unchanged
      C in this process.if the new trust region step turns out to satisfy
      C SNORM > 0.5 * DELTA and VQUAD >= 0 again(i.e., the 'improved'model
      C still suffers from rounding errors), then the program will goto 530
      C and then goto 20, where shifting will not happen because either XBASE
      C was already shifted to the current best point in last step, or XBASE
      C is close to the current best point.Consequently, the model will
      C remain unchanged, and produce the same trust region step again.This
      C leads to an infinite loop.
      C The infinite loop did happen when the MATLAB interface was applied to
      C min atan(x + 100) s.t.x <= -99(x0=-99, npt=3, rhobeg=1, rhoend=1E-6) .
      C The problem does not exist in NEWUOA or BOBYQA, where the program will
      C exit immediately when VQUAD >= 0.
      C To prevent such a loop, here we use IMPRV to record whether the path
      C 530 - - > 20 has already happened for last trust region step.IMPRV = 1
      C implies that last trust region step satisfies VQUAD >= 0 and followed
      C 530 - - > 20.WITH IMPRV = 1, if VQUAD is again nonnegative for the new trust
      C region step, we should not goto 530 but goto 560, where IMPRV will be
      C set to 0 and DELTA will be reduced.Otherwise, an infinite loop would happen.
      C if(KSAVE == 0 .and. VQUAD >= ZERO) goto 530
      if (KSAVE == 0 .and. .not. (VQUAD < ZERO)) then
          if (IMPRV == 1) then
              goto 560
          else
              IMPRV = 1
              goto 530
          end if
      else
          IMPRV = 0
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      C
      C Calculate the next value of the objective function.the difference
      C between the actual new value of F and the value predicted by the
      C model is recorded in DIFF.
      C
220   NF = NF + 1
      if (NF > MAXFUN) then
          NF = NF - 1
          if (IPRINT > 0) print 230
230       format(/4X, 'Return from LINCOA because CALFUN has been',
1         ' called MAXFUN times.')
          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO = 3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          goto 600
      end if
      XDIFF = ZERO
      do I = 1, N
          XNEW(I) = XOPT(I) + STEP(I)
          X(I) = XBASE(I) + XNEW(I)
          XDIFF = XDIFF + (X(I) - XSAV(I))**2
      end do
      XDIFF = DSQRT(XDIFF)
      if (KSAVE == -1) XDIFF = RHO
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C if(XDIFF <= TENTH * RHO .or. XDIFF >= DELTA + DELTA) then
      if (.not. (XDIFF > TENTH * RHO .and. XDIFF < DELTA + DELTA)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IFEAS = 0
          if (IPRINT > 0) print 250
250       format(/4X, 'Return from LINCOA because rounding errors',
1         ' prevent reasonable changes to X.')
          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO = 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          goto 600
      end if
      if (KSAVE <= 0) IFEAS = 1
      F = DFLOAT(IFEAS)
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do I = 1, N
          if (X(I) /= X(I)) then
              F = X(I)  ! Set F to NaN
              if (NF == 1) then
                  FOPT = F
                  XOPT(1:N) = ZERO
              end if
              INFO = -1
              goto 600
          end if
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call CALFUN(N, X, F)
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C By Tom(on 04 - 06 - 2019):
      if (F /= F .or. F > ALMOST_INFINITY) then
          if (NF == 1) then
              FOPT = F
              XOPT(1:N) = ZERO
          end if
          INFO = -2
          goto 600
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (IPRINT == 3) then
          print 260, NF, F, (X(I), I=1, N)
260       format(/4X, 'Function number', I6, '    F =', 1PD18.10,
1         '    The corresponding X is:'/(2X, 5D15.6))
      end if
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C if(KSAVE == -1) goto 600
      if (KSAVE == -1) then
          INFO = 0
          goto 600
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIFF = F - FOPT - VQUAD
      C
      C if X is feasible, then set DFFALT to the difference between the new
      C value of F and the value predicted by the alternative model.
      C
      if (IFEAS == 1 .and. ITEST < 3) then
          do K = 1, NPT
              PQW(K) = ZERO
              W(K) = FVAL(K) - FVAL(KOPT)
          end do
          do J = 1, NPTM
              SUM = ZERO
              do I = 1, NPT
                  SUM = SUM + W(I) * ZMAT(I, J)
              end do
              if (J < IDZ) SUM = -SUM
              do K = 1, NPT
                  PQW(K) = PQW(K) + SUM * ZMAT(K, J)
              end do
          end do
          VQALT = ZERO
          do K = 1, NPT
              SUM = ZERO
              do J = 1, N
                  SUM = SUM + BMAT(K, J) * STEP(J)
              end do
              VQALT = VQALT + SUM * W(K)
              VQALT = VQALT + PQW(K) * SP(NPT + K) * (HALF * SP(NPT + K) + SP(K))
          end do
          DFFALT = F - FOPT - VQALT
      end if
      if (ITEST == 3) then
          DFFALT = DIFF
          ITEST = 0
      end if
      C
      C Pick the next value of DELTA after a trust region step.
      C
      if (KSAVE == 0) then
          RATIO = (F - FOPT) / VQUAD
          if (RATIO <= TENTH) then
              DELTA = HALF * DELTA
          else if (RATIO <= 0.7D0) then
              DELTA = DMAX1(HALF * DELTA, SNORM)
          else
              TEMP = DSQRT(2.0D0) * DELTA
              DELTA = DMAX1(HALF * DELTA, SNORM + SNORM)
              DELTA = DMIN1(DELTA, TEMP)
          end if
          if (DELTA <= 1.4D0 * RHO) DELTA = RHO
      end if
      C
      C Update BMAT, ZMAT and IDZ, so that the KNEW - th interpolation point
      C can be moved.if STEP is a trust region step, then KNEW is zero at
      C present, but a positive value is picked by subroutine UPDATE.
      C
      call UPDATE(N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM, SP, STEP, KOPT,
1     KNEW, PQW, W)
      if (KNEW == 0) then
          if (IPRINT > 0) print 320
320       format(/4X, 'Return from LINCOA because the denominator'
1         ' of the updating formula is zero.')
          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO = 9
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          goto 600
      end if

      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C Zaikun 19 - 03 - 2020
      C exit if BMAT or ZMAT contians NaN
      do J = 1, N
          do I = 1, NDIM
              if (BMAT(I, J) /= BMAT(I, J)) then
                  INFO = -3
                  goto 600
              end if
          end do
      end do
      do J = 1, NPTM
          do I = 1, NPT
              if (ZMAT(I, J) /= ZMAT(I, J)) then
                  INFO = -3
                  goto 600
              end if
          end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      C
      C if ITEST is increased to 3, then the next quadratic model is the
      C one whose second derivative matrix is least subject to the new
      C interpolation conditions.Otherwise the new model is constructed
      C by the symmetric Broyden method in the usual way.
      C
      if (IFEAS == 1) then
          ITEST = ITEST + 1
          if (DABS(DFFALT) >= TENTH * DABS(DIFF)) ITEST = 0
      end if
      C
      C Update the second derivatives of the model by the symmetric Broyden
      C method, using PQW for the second derivative parameters of the new
      C KNEW - th Lagrange function.the contribution from the old parameter
      C PQ(KNEW) is included in the second derivative matrix HQ.W is used
      C later for the gradient of the new KNEW - th Lagrange function.
      C
      if (ITEST < 3) then
          do K = 1, NPT
              PQW(K) = ZERO
          end do
          do J = 1, NPTM
              TEMP = ZMAT(KNEW, J)
              if (TEMP /= ZERO) then
                  if (J < IDZ) TEMP = -TEMP
                  do K = 1, NPT
                      PQW(K) = PQW(K) + TEMP * ZMAT(K, J)
                  end do
              end if
          end do
          IH = 0
          do I = 1, N
              W(I) = BMAT(KNEW, I)
              TEMP = PQ(KNEW) * XPT(KNEW, I)
              do J = 1, I
                  IH = IH + 1
                  HQ(IH) = HQ(IH) + TEMP * XPT(KNEW, J)
              end do
          end do
          PQ(KNEW) = ZERO
          do K = 1, NPT
              PQ(K) = PQ(K) + DIFF * PQW(K)
          end do
      end if
      C
      C include the new interpolation point with the corresponding updates of
      C SP.Also make the changes of the symmetric Broyden method to GOPT at
      C the old XOPT if ITEST is less than 3.
      C
      FVAL(KNEW) = F
      SP(KNEW) = SP(KOPT) + SP(NPT + KOPT)
      SSQ = ZERO
      do I = 1, N
          XPT(KNEW, I) = XNEW(I)
          SSQ = SSQ + STEP(I)**2
      end do
      SP(NPT + KNEW) = SP(NPT + KOPT) + SSQ
      if (ITEST < 3) then
          do K = 1, NPT
              TEMP = PQW(K) * SP(K)
              do I = 1, N
                  W(I) = W(I) + TEMP * XPT(K, I)
              end do
          end do
          do I = 1, N
              GOPT(I) = GOPT(I) + DIFF * W(I)
          end do
      end if
      C
      C Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
      C least calculated value so far with a feasible vector of variables.
      C
      if (F < FOPT .and. IFEAS == 1) then
          FOPT = F
          do J = 1, N
              XSAV(J) = X(J)
              XOPT(J) = XNEW(J)
          end do
          KOPT = KNEW
          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          C By Tom(on 04 - 06 - 2019):
          if (FOPT <= FTARGET) then
              INFO = 1
              goto 616
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SNORM = DSQRT(SSQ)
          do J = 1, M
              if (RESCON(J) >= DELTA + SNORM) then
                  RESCON(J) = SNORM - RESCON(J)
              else
                  RESCON(J) = RESCON(J) + SNORM
                  if (RESCON(J) + DELTA > ZERO) then
                      TEMP = B(J)
                      do I = 1, N
                          TEMP = TEMP - XOPT(I) * AMAT(I, J)
                      end do
                      TEMP = DMAX1(TEMP, ZERO)
                      if (TEMP >= DELTA) TEMP = -TEMP
                      RESCON(J) = TEMP
                  end if
              end if
          end do
          do K = 1, NPT
              SP(K) = SP(K) + SP(NPT + K)
          end do
          C
          C Also revise GOPT when symmetric Broyden updating is applied.
          C
          if (ITEST < 3) then
              IH = 0
              do J = 1, N
                  do I = 1, J
                      IH = IH + 1
                      if (I < J) GOPT(J) = GOPT(J) + HQ(IH) * STEP(I)
                      GOPT(I) = GOPT(I) + HQ(IH) * STEP(J)
                  end do
              end do
              do K = 1, NPT
                  TEMP = PQ(K) * SP(NPT + K)
                  do I = 1, N
                      GOPT(I) = GOPT(I) + TEMP * XPT(K, I)
                  end do
              end do
          end if
      end if
      C
      C Replace the current model by the least Frobenius norm interpolant if
      C this interpolant gives substantial reductions in the predictions
      C of values of F at feasible points.
      C
      if (ITEST == 3) then
          do K = 1, NPT
              PQ(K) = ZERO
              W(K) = FVAL(K) - FVAL(KOPT)
          end do
          do J = 1, NPTM
              SUM = ZERO
              do I = 1, NPT
                  SUM = SUM + W(I) * ZMAT(I, J)
              end do
              if (J < IDZ) SUM = -SUM
              do K = 1, NPT
                  PQ(K) = PQ(K) + SUM * ZMAT(K, J)
              end do
          end do
          do J = 1, N
              GOPT(J) = ZERO
              do I = 1, NPT
                  GOPT(J) = GOPT(J) + W(I) * BMAT(I, J)
              end do
          end do
          do K = 1, NPT
              TEMP = PQ(K) * SP(K)
              do I = 1, N
                  GOPT(I) = GOPT(I) + TEMP * XPT(K, I)
              end do
          end do
          do IH = 1, NH
              HQ(IH) = ZERO
          end do
      end if
      C
      C if a trust region step has provided a sufficient decrease in F, then
      C branch for another trust region calculation.Every iteration that
      C takes a model step is followed by an attempt to take a trust region
      C step.
      C
      KNEW = 0
      if (KSAVE > 0) goto 20
      if (RATIO >= TENTH) goto 20
      C
      C Alternatively, find out if the interpolation points are close enough
      C to the best point so far.
      C
530   DISTSQ = DMAX1(DELTA * DELTA, 4.0D0 * RHO * RHO)
      do K = 1, NPT
          SUM = ZERO
          do J = 1, N
              SUM = SUM + (XPT(K, J) - XOPT(J))**2
          end do
          if (SUM > DISTSQ) then
              KNEW = K
              DISTSQ = SUM
          end if
      end do
      C
      C if KNEW is positive, then branch back for the next iteration, which
      C will generate a "model step".Otherwise, if the current iteration
      C has reduced F, or if DELTA was above its lower bound when the last
      C trust region step was calculated, then try a "trust region"step
      C instead.
      C
      if (KNEW > 0) goto 20
      KNEW = 0
      if (FOPT < FSAVE) goto 20
      if (DELSAV > RHO) goto 20
      C
      C The calculations with the current value of RHO are complete.
      C Pick the next value of RHO.
      C
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C Zaikun 15 - 08 - 2019
      C See the comments below line number 210
      C 560 if(RHO > RHOEND) then
560   IMPRV = 0
      if (RHO > RHOEND) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          DELTA = HALF * RHO
          if (RHO > 250.0D0 * RHOEND) then
              RHO = TENTH * RHO
          else if (RHO <= 16.0D0 * RHOEND) then
              RHO = RHOEND
          else
              RHO = DSQRT(RHO * RHOEND)
          end if
          DELTA = DMAX1(DELTA, RHO)
          if (IPRINT >= 2) then
              if (IPRINT >= 3) print 570
570           format(5X)
              print 580, RHO, NF
580           format(/4X, 'New RHO =', 1PD11.4, 5X, 'Number of',
1             ' function values =', I6)
              print 590, FOPT, (XBASE(I) + XOPT(I), I=1, N)
590           format(4X, 'Least value of F =', 1PD23.15, 9X,
1             'The corresponding X is:'/(2X, 5D15.6))
          end if
          goto 10
      end if
      C
      C return from the calculation, after branching to label 220 for another
      C Newton - Raphson step if it has not been tried before.
      C
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INFO = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (KSAVE == -1) goto 220
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C 600 if(FOPT <= F .or. IFEAS == 0) then
600   if (FOPT <= F .or. IFEAS == 0 .or. F /= F) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              X(I) = XSAV(I)
          end do
          F = FOPT
      end if
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      C if(IPRINT >= 1) then
616   if (IPRINT >= 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          print 620, NF
620       format(/4X, 'At the return from LINCOA', 5X,
1         'Number of function values =', I6)
          print 590, F, (X(I), I=1, N)
      end if
      W(1) = F
      W(2) = DFLOAT(NF) + HALF
      return
      end subroutine LINCOB
