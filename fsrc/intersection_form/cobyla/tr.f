!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of tr.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 07-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      module trustregion_mod

      contains

      subroutine trstlp(n, m, A, b, rho, d, ifull, iact)

! Generic modules
      use consts_mod, only : RP, IK, ZERO, TWO, HALF, TENTH, HUGENUM, DE&
     &BUGGING, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf, is_finite
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : inprod, matprod, eye, planerot, isminor

      implicit none

      integer(IK), intent(in) :: n
      integer(IK), intent(in) :: m
      real(RP), intent(in) :: A(:, :)
      !(n, m+1)
      real(RP), intent(in) :: b(:)
      real(RP), intent(in) :: rho
      real(RP), intent(inout) :: d(:)
      integer(IK), intent(out) :: ifull
      integer(IK), intent(out) :: iact(:)


      real(RP) :: z(n, n)
      real(RP) :: zdota(n)
      real(RP) :: vmultc(m + 1)
      real(RP) :: sdirn(n)
      real(RP) :: vmultd(m + 1)
      real(RP) :: cnew(n)
      real(RP) :: dnew(n)


      real(RP) :: alpha
      real(RP) :: beta
      real(RP) :: dd
      real(RP) :: optnew
      real(RP) :: optold
      real(RP) :: ratio
      real(RP) :: cstrv
      real(RP) :: resold
      real(RP) :: sd
      real(RP) :: sp
      real(RP) :: spabs
      real(RP) :: ss
      real(RP) :: step
      real(RP) :: stpful
      real(RP) :: summ
      real(RP) :: summabs
      real(RP) :: summd
      real(RP) :: temp
      real(RP) :: tempa
      real(RP) :: tot
      real(RP) :: vsave
      real(RP) :: zdotv
      real(RP) :: zdotw
      real(RP) :: zdvabs
      real(RP) :: zdwabs
      real(RP) :: dsav(size(d))
      ! N
      integer(IK) :: i
      integer(IK) :: icon
      integer(IK) :: icount
      integer(IK) :: isave
      integer(IK) :: iterc
      integer(IK) :: j
      integer(IK) :: k
      integer(IK) :: kk
      integer(IK) :: kl
      integer(IK) :: kp
      integer(IK) :: kw
      integer(IK) :: mcon
      integer(IK) :: nact
      integer(IK) :: nactx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine calculates an N-component vector D by the following two stages. In the first
! stage, D is set to the shortest vector that minimizes the greatest violation of the constraints
!       dot_product(A(1:N, K), D) >= B(K),  K= 1, 2, 3, ..., M,
! subject to the Euclidean length of D being at most RHO. If its length is strictly less than RHO,
! then we use the resultant freedom in D to minimize the objective function
!       dot_product(-A(1:N, M+1), D)
! subject to no increase in any greatest constraint violation. This notation allows the gradient of
! the objective function to be regarded as the gradient of a constraint. Therefore the two stages
! are distinguished by MCON == M and MCON > M respectively.
!
! It is possible that a degeneracy may prevent D from attaining the target length RHO. Then the
! value IFULL=0 would be set, but usually IFULL=1 on return.
!
! In general NACT is the number of constraints in the active set and IACT(1),...,IACT(NACT) are
! their indices, while the remainder of IACT contains a permutation of the remaining constraint
! indices.
!
! Further, Z is an orthogonal matrix whose first NACT columns can be regarded as the result of
! Gram-Schmidt applied to the active constraint gradients. For J=1, 2, ..., NACT, the number
! ZDOTA(J) is the scalar product of the J-th column of Z with the gradient of the J-th active
! constraint. D is the current vector of variables and here the residuals of the active constraints
! should be zero. Further, the active constraints have nonnegative Lagrange multipliers that are
! held at the beginning of VMULTC. The remainder of this vector holds the residuals of the inactive
! constraints at d, the ordering of the components of VMULTC being in agreement with the permutation
! of the indices of the constraints that is in IACT. All these residuals are nonnegative, which is
! achieved by the shift CSTRV that makes the least residual zero.
!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
      ITERC = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IFULL = 1
      MCON = M
      NACT = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      CSTRV=0.0
      CSTRV = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! Question: what is the size of B? M, M+1, or M+2?
! It seems to be M+1. Then when is B(M+1) used?
! What is the size of iact, vmultc? m or m+1?

! Initialize Z and some other variables. The value of CSTRV will be appropriate to D=0, while ICON
! will be the index of a most violated constraint if CSTRV is positive. Usually during the first
! stage the vector SDIRN gives a search direction that reduces all the active constraint violations
! by one simultaneously.
      z = eye(n, n)
      d = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      cstrv = maxval([b(1:m), ZERO])
      icon = maxloc(b(1:m), dim=1)
      iact(1:m) = [(k, k=1, m)]
      vmultc(1:m) = cstrv - b(1:m)

      if (cstrv <= zero) goto 480

      sdirn = ZERO

      dsav = d
      !!!! HOW TO AVOID THIS???


! End the current stage of the calculation if 3 consecutive iterations have either failed to reduce
! the best calculated value of the objective function or to increase the number of active
! constraints since the best value was calculated. This strategy prevents cycling, but there is
! a remote possibility that it will cause premature termination.
!
60    optold = ZERO

      icount = 0_IK
70    if (mcon == m) then
          optnew = cstrv
      else
          optnew = -inprod(d, A(:, mcon))
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019
! The original code can still encounter infinite cycling,
! which did happen when testing the following CUTEst problems:
! DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely
! large values in A (up to 10^219).
! To avoid wasting energy, we do the following.

      if (is_finite(sum(abs(d)))) then
          dsav = d
      else
          d = dsav
          return
      end if
      iterc = iterc + 1
      if (iterc > min(10000_IK, 100_IK * n)) then
          return
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! What is ICOUNT ????
      if (icount == 0 .or. optnew < optold) then
          optold = optnew
          nactx = nact
          icount = 3
      else if (nact > nactx) then
          nactx = nact
          icount = 3
      else
          icount = icount - 1
          if (icount == 0) goto 490
      end if


! If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to the active set. Apply
! Givens rotations so that the last N-NACT-1 columns of Z are orthogonal to the gradient of the new
! constraint, a scalar product being set to zero if its nonzero value could be due to computer
! rounding errors.
      if (icon <= nact) goto 260

      kk = iact(icon)
      cnew = A(:, kk)
      !!! KK = ???

      tot = ZERO
      do k = n, nact + 1, -1
          sp = inprod(z(:, k), cnew)
          if (isminor(sp, sum(abs(z(:, k) * cnew)))) then
              sp = ZERO
          end if
          if (abs(tot) <= ZERO) then
              tot = sp
          else
              kp = k + 1_IK
              z(:, [k, kp]) = matprod(z(:, [k, kp]), planerot([sp, tot])&
     &)
              tot = sqrt(sp * sp + tot * tot)
! After the rotation,
! Z_NEW(:, K)'*CNEW = SQRT(Z(:, K)'*CNEW^2 + Z(:, KP)'*CNEW^2),
! Z_NEW(:, KP)'*CNEW = 0
          end if
      end do
!
!     Add the new constraint if this can be done without a deletion from the
!     active set.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TOT .NE. 0.0) THEN
      if (TOT /= 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          NACT = NACT + 1
          ZDOTA(NACT) = TOT
          VMULTC(ICON) = VMULTC(NACT)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTC(NACT)=0.0
          VMULTC(NACT) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          goto 210
      end if
!
!     The next instruction is reached if a deletion has to be made from the
!     active set in order to make room for the new active constraint, because
!     the new constraint gradient is a linear combination of the gradients of
!     the old active constraints. Set the elements of VMULTD to the multipliers
!     of the linear combination. Further, set IOUT to the index of the
!     constraint to be deleted, but branch if no suitable index can be found.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=-1.0
      RATIO = -1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      K = NACT
      !!!! What if K = 0?
! It seems to lead to a memory error when VMULTD(K) is accessed (by why not Z(I, K) and others?????)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  130 ZDOTV=0.0
!      ZDVABS=0.0
130   ZDOTV = 0.0D0
      ZDVABS = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          TEMP = Z(I, K) * cnew(I)
          ZDOTV = ZDOTV + TEMP
          ZDVABS = ZDVABS + DABS(TEMP)
      end do
      if (.not. isminor(zdotv, zdvabs)) then
          TEMP = ZDOTV / ZDOTA(K)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .GT. 0.0 .AND. IACT(K) .LE. M) THEN
          if (TEMP > 0.0D0 .and. IACT(K) <= M) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TEMPA = VMULTC(K) / TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              IF (RATIO .LT. 0.0 .OR. TEMPA .LT. RATIO) THEN
              if (RATIO < 0.0D0 .or. TEMPA < RATIO) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  RATIO = TEMPA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: IOUT is never used
!                  IOUT=K
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              end if
          end if
          if (K >= 2) then
              KW = IACT(K)
              do I = 1, N
                  cnew(I) = cnew(I) - TEMP * A(I, KW)
              end do
          end if
          VMULTD(K) = TEMP
      else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTD(K)=0.0
!!!!! Zaikun 2021-06-27 XXXXXXXXXXXXXXXXXX
!    VMULTD(K) = 0.0D0  !!!! This seems to lead to a memory error ("free() invalid pointer") when
!    k = 0.
          if (K > 0) VMULTD(K) = 0.0D0
!!!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
      K = K - 1
      if (K > 0) goto 130
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (RATIO .LT. 0.0) GOTO 490
      if (RATIO < 0.0D0) goto 490
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Revise the Lagrange multipliers and reorder the active constraints so
!     that the one to be replaced is at the end of the list. Also calculate the
!     new value of ZDOTA(NACT) and branch if it is not acceptable.
!
      do K = 1, NACT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  160 VMULTC(K)=AMAX1(0.0,VMULTC(K)-RATIO*VMULTD(K))
          VMULTC(K) = DMAX1(0.0D0, VMULTC(K) - RATIO * VMULTD(K))
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ICON < NACT) then
          ISAVE = IACT(ICON)
          VSAVE = VMULTC(ICON)
          K = ICON
170   KP = K + 1
          KW = IACT(KP)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
          SP = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              SP = SP + Z(I, K) * A(I, KW)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
          TEMP = DSQRT(SP * SP + ZDOTA(KP)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ALPHA = ZDOTA(KP) / TEMP
          BETA = SP / TEMP
          ZDOTA(KP) = ALPHA * ZDOTA(K)
          ZDOTA(K) = TEMP
          do I = 1, N
              TEMP = ALPHA * Z(I, KP) + BETA * Z(I, K)
              Z(I, KP) = ALPHA * Z(I, K) - BETA * Z(I, KP)
              Z(I, K) = TEMP
          end do
          IACT(K) = KW
          VMULTC(K) = VMULTC(KP)
          K = KP
          if (K < NACT) goto 170
          IACT(K) = ISAVE
          VMULTC(K) = VSAVE
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
      TEMP = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          TEMP = TEMP + Z(I, NACT) * A(I, KK)
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TEMP .EQ. 0.0) GOTO 490
      if (TEMP == 0.0D0) goto 490
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ZDOTA(NACT) = TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VMULTC(ICON)=0.0
      VMULTC(ICON) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VMULTC(NACT) = RATIO
!
!     Update IACT and ensure that the objective function continues to be
!     treated as the last active constraint when MCON>M.
!
210   IACT(ICON) = IACT(NACT)
      IACT(NACT) = KK
      if (MCON > M .and. KK /= MCON) then
          K = NACT - 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
          SP = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              SP = SP + Z(I, K) * A(I, KK)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(NACT)**2)
          TEMP = DSQRT(SP * SP + ZDOTA(NACT)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ALPHA = ZDOTA(NACT) / TEMP
          BETA = SP / TEMP
          ZDOTA(NACT) = ALPHA * ZDOTA(K)
          ZDOTA(K) = TEMP
          do I = 1, N
              TEMP = ALPHA * Z(I, NACT) + BETA * Z(I, K)
              Z(I, NACT) = ALPHA * Z(I, K) - BETA * Z(I, NACT)
              Z(I, K) = TEMP
          end do
          IACT(NACT) = IACT(K)
          IACT(K) = KK
          TEMP = VMULTC(K)
          VMULTC(K) = VMULTC(NACT)
          VMULTC(NACT) = TEMP
      end if
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
      if (MCON > M) goto 320
      KK = IACT(NACT)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
      TEMP = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          TEMP = TEMP + SDIRN(I) * A(I, KK)
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=TEMP-1.0
      TEMP = TEMP - 1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TEMP = TEMP / ZDOTA(NACT)
      do I = 1, N
          SDIRN(I) = SDIRN(I) - TEMP * Z(I, NACT)
      end do
      goto 340
!
!     Delete the constraint that has the index IACT(ICON) from the active set.
!
260   if (ICON < NACT) then
          ISAVE = IACT(ICON)
          VSAVE = VMULTC(ICON)
          K = ICON
270   KP = K + 1
          KK = IACT(KP)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
          SP = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              SP = SP + Z(I, K) * A(I, KK)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
          TEMP = DSQRT(SP * SP + ZDOTA(KP)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ALPHA = ZDOTA(KP) / TEMP
          BETA = SP / TEMP
          ZDOTA(KP) = ALPHA * ZDOTA(K)
          ZDOTA(K) = TEMP
          do I = 1, N
              TEMP = ALPHA * Z(I, KP) + BETA * Z(I, K)
              Z(I, KP) = ALPHA * Z(I, K) - BETA * Z(I, KP)
              Z(I, K) = TEMP
          end do
          IACT(K) = KK
          VMULTC(K) = VMULTC(KP)
          K = KP
          if (K < NACT) goto 270
          IACT(K) = ISAVE
          VMULTC(K) = VSAVE
      end if
      NACT = NACT - 1
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
      if (MCON > M) goto 320
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
      TEMP = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          TEMP = TEMP + SDIRN(I) * Z(I, NACT + 1)
      end do
      do I = 1, N
          SDIRN(I) = SDIRN(I) - TEMP * Z(I, NACT + 1)
      end do
      GO TO 340
!
!     Pick the next search direction of stage two.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  320 TEMP=1.0/ZDOTA(NACT)
320   TEMP = 1.0D0 / ZDOTA(NACT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          SDIRN(I) = TEMP * Z(I, NACT)
      end do
!
!     Calculate the step to the boundary of the trust region or take the step
!     that reduces CSTRV to zero. The two statements below that include the
!     factor 1.0E-6 prevent some harmless underflows that occurred in a test
!     calculation. Further, we skip the step if it could be zero within a
!     reasonable tolerance for computer rounding errors.
!
340   DD = RHO * RHO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SD=0.0
!      SS=0.0
      SD = 0.0D0
      SS = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (ABS(d(I)) .GE. 1.0E-6*RHO) DD=DD-d(I)**2
          if (DABS(d(I)) >= 1.0D-6 * RHO) DD = DD - d(I)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SD = SD + d(I) * SDIRN(I)
          SS = SS + SDIRN(I)**2
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (DD .LE. 0.0) GOTO 490
!      TEMP=SQRT(SS*DD)
!      IF (ABS(SD) .GE. 1.0E-6*TEMP) TEMP=SQRT(SS*DD+SD*SD)
      if (DD <= 0.0D0) goto 490
      TEMP = DSQRT(SS * DD)
      if (DABS(SD) >= 1.0D-6 * TEMP) TEMP = DSQRT(SS * DD + SD * SD)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      STPFUL = DD / (TEMP + SD)
      STEP = STPFUL
      if (MCON == M) then
          if (isminor(cstrv, step)) then
              goto 480
          end if
          STEP = DMIN1(STEP, CSTRV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
!
!     Set dnew to the new variables if STEP is the steplength, and reduce
!     CSTRV to the corresponding maximum residual if stage one is being done.
!     Because dnew will be changed during the calculation of some Lagrange
!     multipliers, it will be restored to the following value later.
!
      do I = 1, N
          dnew(I) = d(I) + STEP * SDIRN(I)
      end do
      if (MCON == M) then
          RESOLD = CSTRV
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          CSTRV=0.0
          CSTRV = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do K = 1, NACT
              KK = IACT(K)
              TEMP = B(KK)
              do I = 1, N
                  TEMP = TEMP - A(I, KK) * dnew(I)
              end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          CSTRV=AMAX1(CSTRV,TEMP)
              CSTRV = DMAX1(CSTRV, TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          end do
      end if
!
!     Set VMULTD to the VMULTC vector that would occur if d became dnew. A
!     device is included to force VMULTD(K)=0.0 if deviations from this value
!     can be attributed to computer rounding errors. First calculate the new
!     Lagrange multipliers.
!
      K = NACT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  390 ZDOTW=0.0
!      ZDWABS=0.0
390   ZDOTW = 0.0D0
      ZDWABS = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          TEMP = Z(I, K) * dnew(I)
          ZDOTW = ZDOTW + TEMP
          ZDWABS = ZDWABS + DABS(TEMP)
      end do
      if (isminor(zdotw, zdwabs)) then
          zdotw = ZERO
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VMULTD(K) = ZDOTW / ZDOTA(K)
      if (K >= 2) then
          KK = IACT(K)
          do I = 1, N
              dnew(I) = dnew(I) - VMULTD(K) * A(I, KK)
          end do
          K = K - 1
          goto 390
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (MCON .GT. M) VMULTD(NACT)=AMAX1(0.0,VMULTD(NACT))
      if (MCON > M) VMULTD(NACT) = DMAX1(0.0D0, VMULTD(NACT))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Complete VMULTC by finding the new constraint residuals.
!
      do I = 1, N
          dnew(I) = d(I) + STEP * SDIRN(I)
      end do
      if (MCON > NACT) then
          KL = NACT + 1
          do K = KL, MCON
              KK = IACT(K)
              SUMM = CSTRV - B(KK)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SUMMABS=CSTRV+ABS(B(KK))
              SUMMABS = CSTRV + DABS(B(KK))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do I = 1, N
                  TEMP = A(I, KK) * dnew(I)
                  SUMM = SUMM + TEMP
                  SUMMABS = SUMMABS + DABS(TEMP)
              end do
              if (isminor(summ, summabs)) then
                  summ = ZERO
              end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              VMULTD(K) = SUMM
          end do
      end if
!
!     Calculate the fraction of the step from d to dnew that will be taken.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=1.0
      RATIO = 1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICON = 0
      do K = 1, MCON
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (VMULTD(K) .LT. 0.0) THEN
          if (VMULTD(K) < 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TEMP = VMULTC(K) / (VMULTC(K) - VMULTD(K))
              if (TEMP < RATIO) then
                  RATIO = TEMP
                  ICON = K
              end if
          end if
      end do
!
!     Update d, VMULTC and CSTRV.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=1.0-RATIO
      TEMP = 1.0D0 - RATIO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          d(I) = TEMP * d(I) + RATIO * dnew(I)
      end do
      do K = 1, MCON
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  470 VMULTC(K)=AMAX1(0.0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
          VMULTC(K) = DMAX1(0.0D0, TEMP * VMULTC(K) + RATIO * VMULTD(K))
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (MCON == M) CSTRV = RESOLD + RATIO * (CSTRV - RESOLD)
!
!     If the full step is not acceptable then begin another iteration.
!     Otherwise switch to stage two or end the calculation.
!
      if (ICON > 0) goto 70
      if (STEP == STPFUL) goto 500
480   MCON = M + 1
      ICON = MCON
      IACT(MCON) = MCON
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VMULTC(MCON)=0.0
      VMULTC(MCON) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      goto 60
!
!     We employ any freedom that may be available to reduce the objective
!     function before returning a d whose length is less than RHO.
!
490   if (MCON == M) goto 480
      IFULL = 0
500   return

      end subroutine trstlp

      end module trustregion_mod