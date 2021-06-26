!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of trustregion.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 26-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      module trustregion_mod

      contains

      subroutine trstlp(n, m, A, B, rho, dx, ifull, iact)

! Generic modules
      use consts_mod, only : RP, IK, ZERO, TWO, HALF, TENTH, HUGENUM, DE&
     &BUGGING, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : calquad, inprod

      implicit none

      integer(IK), intent(IN) :: N
      integer(IK), intent(IN) :: M
      real(RP), intent(IN) :: A(:, :)
      !(N, :)
      real(RP), intent(IN) :: B(:)
      real(RP), intent(IN) :: Rho
      real(RP), intent(INOUT) :: Dx(:)
      integer(IK), intent(OUT) :: Ifull
      integer(IK), intent(INOUT) :: Iact(:)


      real(RP) :: Z(N, N)
      real(RP) :: Zdota(N)
      real(RP) :: Vmultc(M + 1)
      real(RP) :: Sdirn(N)
      real(RP) :: Dxnew(N)
      real(RP) :: Vmultd(M)

      real(RP) :: acca
      real(RP) :: accb
      real(RP) :: alpha
      real(RP) :: beta
      real(RP) :: dd
      real(RP) :: optnew
      real(RP) :: optold
      real(RP) :: ratio
      real(RP) :: resmax
      real(RP) :: resold
      real(RP) :: sd
      real(RP) :: sp
      real(RP) :: spabs
      real(RP) :: ss
      real(RP) :: step
      real(RP) :: stpful
      real(RP) :: sum
      real(RP) :: sumabs
      real(RP) :: sumd
      real(RP) :: temp
      real(RP) :: tempa
      real(RP) :: tot
      real(RP) :: vsave
      real(RP) :: zdotv
      real(RP) :: zdotw
      real(RP) :: zdvabs
      real(RP) :: zdwabs
      real(RP) :: dsav(size(dx))
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
!     This subroutine calculates an N-component vector DX by applying the
!     following two stages. In the first stage, DX is set to the shortest
!     vector that minimizes the greatest violation of the constraints
!       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M,
!     subject to the Euclidean length of DX being at most RHO. If its length is
!     strictly less than RHO, then we use the resultant freedom in DX to
!     minimize the objective function
!              -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N)
!     subject to no increase in any greatest constraint violation. This
!     notation allows the gradient of the objective function to be regarded as
!     the gradient of a constraint. Therefore the two stages are distinguished
!     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
!     degeneracy may prevent DX from attaining the target length RHO. Then the
!     value IFULL=0 would be set, but usually IFULL=1 on return.
!
!     In general NACT is the number of constraints in the active set and
!     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
!     contains a permutation of the remaining constraint indices. Further, Z is
!     an orthogonal matrix whose first NACT columns can be regarded as the
!     result of Gram-Schmidt applied to the active constraint gradients. For
!     J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
!     column of Z with the gradient of the J-th active constraint. DX is the
!     current vector of variables and here the residuals of the active
!     constraints should be zero. Further, the active constraints have
!     nonnegative Lagrange multipliers that are held at the beginning of
!     VMULTC. The remainder of this vector holds the residuals of the inactive
!     constraints at DX, the ordering of the components of VMULTC being in
!     agreement with the permutation of the indices of the constraints that is
!     in IACT. All these residuals are nonnegative, which is achieved by the
!     shift RESMAX that makes the least residual zero.
!
!     Initialize Z and some other variables. The value of RESMAX will be
!     appropriate to DX=0, while ICON will be the index of a most violated
!     constraint if RESMAX is positive. Usually during the first stage the
!     vector SDIRN gives a search direction that reduces all the active
!     constraint violations by one simultaneously.
!
!

      write (10, *) n, m, A, B, rho, dx, ifull, iact

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
      ITERC = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IFULL = 1
      MCON = M
      NACT = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RESMAX=0.0
      RESMAX = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          do J = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   10 Z(I,J)=0.0
!      Z(I,I)=1.0
!   20 DX(I)=0.0
              Z(I, J) = 0.0D0
          end do
          Z(I, I) = 1.0D0
          DX(I) = 0.0D0
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (M >= 1) then
          do K = 1, M
              if (B(K) > RESMAX) then
                  RESMAX = B(K)
                  ICON = K
              end if
          end do
          do K = 1, M
              IACT(K) = K
              VMULTC(K) = RESMAX - B(K)
          end do
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (RESMAX .EQ. 0.0) GOTO 480
      if (RESMAX == 0.0D0) goto 480
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   50 SDIRN(I)=0.0
          SDIRN(I) = 0.0D0
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
      do I = 1, N
          DSAV(I) = DX(I)
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     End the current stage of the calculation if 3 consecutive iterations
!     have either failed to reduce the best calculated value of the objective
!     function or to increase the number of active constraints since the best
!     value was calculated. This strategy prevents cycling, but there is a
!     remote possibility that it will cause premature termination.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   60 OPTOLD=0.0
60    OPTOLD = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICOUNT = 0
70    if (MCON == M) then
          OPTNEW = RESMAX
      else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          OPTNEW=0.0
          OPTNEW = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              OPTNEW = OPTNEW - DX(I) * A(I, MCON)
          end do
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019
! The original code can still encounter infinite cycling,
! which did happen when testing the following CUTEst problems:
! DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322.
! Indeed, in all these cases, Inf and NaN appear in D due to extremely
! large values in A (up to 10^219).
! To avoid wasting energy, we do the following.
      SUMD = 0.0D0
      do I = 1, N
          SUMD = SUMD + DABS(DX(I))
      end do
      if (SUMD >= 1.0D100 .or. SUMD /= SUMD) then
          do I = 1, N
              DX(I) = DSAV(I)
          end do
          goto 500
      else
          do I = 1, N
              DSAV(I) = DX(I)
          end do
      end if
      ITERC = ITERC + 1
      if (ITERC > min(10000, 100 * N)) then
          goto 500
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ICOUNT == 0 .or. OPTNEW < OPTOLD) then
          OPTOLD = OPTNEW
          NACTX = NACT
          ICOUNT = 3
      else if (NACT > NACTX) then
          NACTX = NACT
          ICOUNT = 3
      else
          ICOUNT = ICOUNT - 1
          if (ICOUNT == 0) goto 490
      end if
!
!     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
!     the active set. Apply Givens rotations so that the last N-NACT-1 columns
!     of Z are orthogonal to the gradient of the new constraint, a scalar
!     product being set to zero if its nonzero value could be due to computer
!     rounding errors. The array DXNEW is used for working space.
!
      if (ICON <= NACT) goto 260
      KK = IACT(ICON)
      do I = 1, N
          DXNEW(I) = A(I, KK)
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TOT=0.0
      TOT = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      K = N
100   if (K > NACT) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
!          SPABS=0.0
          SP = 0.0D0
          SPABS = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              TEMP = Z(I, K) * DXNEW(I)
              SP = SP + TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  110     SPABS=SPABS+ABS(TEMP)
!          ACCA=SPABS+0.1*ABS(SP)
!          ACCB=SPABS+0.2*ABS(SP)
!          IF (SPABS .GE. ACCA .OR. ACCA .GE. ACCB) SP=0.0
!          IF (TOT .EQ. 0.0) THEN
              SPABS = SPABS + DABS(TEMP)
          end do
          ACCA = SPABS + 0.1D0 * DABS(SP)
          ACCB = SPABS + 0.2D0 * DABS(SP)
          if (SPABS >= ACCA .or. ACCA >= ACCB) SP = 0.0D0
          if (TOT == 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              TOT = SP
          else
              KP = K + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              TEMP=SQRT(SP*SP+TOT*TOT)
              TEMP = DSQRT(SP * SP + TOT * TOT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ALPHA = SP / TEMP
              BETA = TOT / TEMP
              TOT = TEMP
              do I = 1, N
                  TEMP = ALPHA * Z(I, K) + BETA * Z(I, KP)
                  Z(I, KP) = ALPHA * Z(I, KP) - BETA * Z(I, K)
                  Z(I, K) = TEMP
              end do
          end if
          K = K - 1
          goto 100
      end if
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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  130 ZDOTV=0.0
!      ZDVABS=0.0
130   ZDOTV = 0.0D0
      ZDVABS = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          TEMP = Z(I, K) * DXNEW(I)
          ZDOTV = ZDOTV + TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  140 ZDVABS=ZDVABS+ABS(TEMP)
!      ACCA=ZDVABS+0.1*ABS(ZDOTV)
!      ACCB=ZDVABS+0.2*ABS(ZDOTV)
          ZDVABS = ZDVABS + DABS(TEMP)
      end do
      ACCA = ZDVABS + 0.1D0 * DABS(ZDOTV)
      ACCB = ZDVABS + 0.2D0 * DABS(ZDOTV)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ZDVABS < ACCA .and. ACCA < ACCB) then
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
                  DXNEW(I) = DXNEW(I) - TEMP * A(I, KW)
              end do
          end if
          VMULTD(K) = TEMP
      else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTD(K)=0.0
          VMULTD(K) = 0.0D0
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
!     that reduces RESMAX to zero. The two statements below that include the
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
!      IF (ABS(DX(I)) .GE. 1.0E-6*RHO) DD=DD-DX(I)**2
          if (DABS(DX(I)) >= 1.0D-6 * RHO) DD = DD - DX(I)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SD = SD + DX(I) * SDIRN(I)
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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          ACCA=STEP+0.1*RESMAX
!          ACCB=STEP+0.2*RESMAX
          ACCA = STEP + 0.1D0 * RESMAX
          ACCB = STEP + 0.2D0 * RESMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (STEP >= ACCA .or. ACCA >= ACCB) goto 480
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          STEP=AMIN1(STEP,RESMAX)
          STEP = DMIN1(STEP, RESMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
!
!     Set DXNEW to the new variables if STEP is the steplength, and reduce
!     RESMAX to the corresponding maximum residual if stage one is being done.
!     Because DXNEW will be changed during the calculation of some Lagrange
!     multipliers, it will be restored to the following value later.
!
      do I = 1, N
          DXNEW(I) = DX(I) + STEP * SDIRN(I)
      end do
      if (MCON == M) then
          RESOLD = RESMAX
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RESMAX=0.0
          RESMAX = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do K = 1, NACT
              KK = IACT(K)
              TEMP = B(KK)
              do I = 1, N
                  TEMP = TEMP - A(I, KK) * DXNEW(I)
              end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RESMAX=AMAX1(RESMAX,TEMP)
              RESMAX = DMAX1(RESMAX, TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          end do
      end if
!
!     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
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
          TEMP = Z(I, K) * DXNEW(I)
          ZDOTW = ZDOTW + TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  400 ZDWABS=ZDWABS+ABS(TEMP)
!      ACCA=ZDWABS+0.1*ABS(ZDOTW)
!      ACCB=ZDWABS+0.2*ABS(ZDOTW)
!      IF (ZDWABS .GE. ACCA .OR. ACCA .GE. ACCB) ZDOTW=0.0
          ZDWABS = ZDWABS + DABS(TEMP)
      end do
      ACCA = ZDWABS + 0.1D0 * DABS(ZDOTW)
      ACCB = ZDWABS + 0.2D0 * DABS(ZDOTW)
      if (ZDWABS >= ACCA .or. ACCA >= ACCB) ZDOTW = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      VMULTD(K) = ZDOTW / ZDOTA(K)
      if (K >= 2) then
          KK = IACT(K)
          do I = 1, N
              DXNEW(I) = DXNEW(I) - VMULTD(K) * A(I, KK)
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
          DXNEW(I) = DX(I) + STEP * SDIRN(I)
      end do
      if (MCON > NACT) then
          KL = NACT + 1
          do K = KL, MCON
              KK = IACT(K)
              SUM = RESMAX - B(KK)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SUMABS=RESMAX+ABS(B(KK))
              SUMABS = RESMAX + DABS(B(KK))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do I = 1, N
                  TEMP = A(I, KK) * DXNEW(I)
                  SUM = SUM + TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  430     SUMABS=SUMABS+ABS(TEMP)
!          ACCA=SUMABS+0.1*ABS(SUM)
!          ACCB=SUMABS+0.2*ABS(SUM)
!          IF (SUMABS .GE. ACCA .OR. ACCA .GE. ACCB) SUM=0.0
                  SUMABS = SUMABS + DABS(TEMP)
              end do
              ACCA = SUMABS + 0.1D0 * DABS(SUM)
              ACCB = SUMABS + 0.2D0 * DABS(SUM)
              if (SUMABS >= ACCA .or. ACCA >= ACCB) SUM = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              VMULTD(K) = SUM
          end do
      end if
!
!     Calculate the fraction of the step from DX to DXNEW that will be taken.
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
!     Update DX, VMULTC and RESMAX.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=1.0-RATIO
      TEMP = 1.0D0 - RATIO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          DX(I) = TEMP * DX(I) + RATIO * DXNEW(I)
      end do
      do K = 1, MCON
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  470 VMULTC(K)=AMAX1(0.0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
          VMULTC(K) = DMAX1(0.0D0, TEMP * VMULTC(K) + RATIO * VMULTD(K))
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (MCON == M) RESMAX = RESOLD + RATIO * (RESMAX - RESOLD)
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
!     function before returning a DX whose length is less than RHO.
!
490   if (MCON == M) goto 480
      IFULL = 0
500   write (10, *) 'ot'
      close (10)
      return

      end subroutine trstlp
      end module trustregion_mod