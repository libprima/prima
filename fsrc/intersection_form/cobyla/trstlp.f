!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of trstlp.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 25-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine trstlp(n, m, A, B, rho, dx, ifull, iact, z, zdota, vmul&
     &tc, sdirn, dxnew, vmultd)

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
      real(RP), intent(INOUT) :: Z(:, :)
      ! (N, :)
      real(RP), intent(INOUT) :: Zdota(:)
      real(RP), intent(INOUT) :: Vmultc(:)
      real(RP), intent(INOUT) :: Sdirn(:)
      real(RP), intent(INOUT) :: Dxnew(:)
      real(RP), intent(INOUT) :: Vmultd(:)

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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
      iterc = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Ifull = 1
      mcon = M
      nact = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RESMAX=0.0
      resmax = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, N
          do j = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   10 Z(I,J)=0.0
!      Z(I,I)=1.0
!   20 DX(I)=0.0
              Z(i, j) = 0.0D0
          end do
          Z(i, i) = 1.0D0
          Dx(i) = 0.0D0
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (M >= 1) then
          do k = 1, M
              if (B(k) > resmax) then
                  resmax = B(k)
                  icon = k
              end if
          end do
          do k = 1, M
              Iact(k) = k
              Vmultc(k) = resmax - B(k)
          end do
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (RESMAX .EQ. 0.0) GOTO 480
      if (resmax == 0.0D0) goto 600
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   50 SDIRN(I)=0.0
          Sdirn(i) = 0.0D0
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 26-06-2019: See the code below line number 80
      do i = 1, N
          dsav(i) = Dx(i)
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
100   optold = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      icount = 0
200   if (mcon == M) then
          optnew = resmax
      else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          OPTNEW=0.0
          optnew = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i = 1, N
              optnew = optnew - Dx(i) * A(i, mcon)
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
      sumd = 0.0D0
      do i = 1, N
          sumd = sumd + DABS(Dx(i))
      end do
      if (sumd >= 1.0D100 .or. sumd /= sumd) then
          do i = 1, N
              Dx(i) = dsav(i)
          end do
          goto 99999
      else
          do i = 1, N
              dsav(i) = Dx(i)
          end do
      end if
      iterc = iterc + 1
      if (iterc > min(10000, 100 * N)) goto 99999
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (icount == 0 .or. optnew < optold) then
          optold = optnew
          nactx = nact
          icount = 3
      elseif (nact > nactx) then
          nactx = nact
          icount = 3
      else
          icount = icount - 1
          if (icount == 0) goto 700
      end if
!
!     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
!     the active set. Apply Givens rotations so that the last N-NACT-1 columns
!     of Z are orthogonal to the gradient of the new constraint, a scalar
!     product being set to zero if its nonzero value could be due to computer
!     rounding errors. The array DXNEW is used for working space.
!
      if (icon <= nact) then
!
!     Delete the constraint that has the index IACT(ICON) from the active set.
!
          if (icon < nact) then
              isave = Iact(icon)
              vsave = Vmultc(icon)
              k = icon
              do
                  kp = k + 1
                  kk = Iact(kp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
                  sp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do i = 1, N
                      sp = sp + Z(i, k) * A(i, kk)
                  end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
                  temp = DSQRT(sp * sp + Zdota(kp)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  alpha = Zdota(kp) / temp
                  beta = sp / temp
                  Zdota(kp) = alpha * Zdota(k)
                  Zdota(k) = temp
                  do i = 1, N
                      temp = alpha * Z(i, kp) + beta * Z(i, k)
                      Z(i, kp) = alpha * Z(i, k) - beta * Z(i, kp)
                      Z(i, k) = temp
                  end do
                  Iact(k) = kk
                  Vmultc(k) = Vmultc(kp)
                  k = kp
                  if (k >= nact) then
                      Iact(k) = isave
                      Vmultc(k) = vsave
                      exit
                  end if
              end do
          end if
          nact = nact - 1
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
          if (mcon > M) goto 400
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
          temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i = 1, N
              temp = temp + Sdirn(i) * Z(i, nact + 1)
          end do
          do i = 1, N
              Sdirn(i) = Sdirn(i) - temp * Z(i, nact + 1)
          end do
          goto 500
      else
          kk = Iact(icon)
          do i = 1, N
              Dxnew(i) = A(i, kk)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TOT=0.0
          tot = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          k = N
          do
              if (k > nact) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
!          SPABS=0.0
                  sp = 0.0D0
                  spabs = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do i = 1, N
                      temp = Z(i, k) * Dxnew(i)
                      sp = sp + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  110     SPABS=SPABS+ABS(TEMP)
!          ACCA=SPABS+0.1*ABS(SP)
!          ACCB=SPABS+0.2*ABS(SP)
!          IF (SPABS .GE. ACCA .OR. ACCA .GE. ACCB) SP=0.0
!          IF (TOT .EQ. 0.0) THEN
                      spabs = spabs + DABS(temp)
                  end do
                  acca = spabs + 0.1D0 * DABS(sp)
                  accb = spabs + 0.2D0 * DABS(sp)
                  if (spabs >= acca .or. acca >= accb) sp = 0.0D0
                  if (tot == 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      tot = sp
                  else
                      kp = k + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              TEMP=SQRT(SP*SP+TOT*TOT)
                      temp = DSQRT(sp * sp + tot * tot)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      alpha = sp / temp
                      beta = tot / temp
                      tot = temp
                      do i = 1, N
                          temp = alpha * Z(i, k) + beta * Z(i, kp)
                          Z(i, kp) = alpha * Z(i, kp) - beta * Z(i, k)
                          Z(i, k) = temp
                      end do
                  end if
                  k = k - 1
                  cycle
              end if
!
!     Add the new constraint if this can be done without a deletion from the
!     active set.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TOT .NE. 0.0) THEN
              if (tot /= 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  nact = nact + 1
                  Zdota(nact) = tot
                  Vmultc(icon) = Vmultc(nact)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTC(NACT)=0.0
                  Vmultc(nact) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  goto 300
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
              ratio = -1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              k = nact
              exit
          end do
          do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  130 ZDOTV=0.0
!      ZDVABS=0.0
              zdotv = 0.0D0
              zdvabs = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do i = 1, N
                  temp = Z(i, k) * Dxnew(i)
                  zdotv = zdotv + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  140 ZDVABS=ZDVABS+ABS(TEMP)
!      ACCA=ZDVABS+0.1*ABS(ZDOTV)
!      ACCB=ZDVABS+0.2*ABS(ZDOTV)
                  zdvabs = zdvabs + DABS(temp)
              end do
              acca = zdvabs + 0.1D0 * DABS(zdotv)
              accb = zdvabs + 0.2D0 * DABS(zdotv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (zdvabs < acca .and. acca < accb) then
                  temp = zdotv / Zdota(k)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .GT. 0.0 .AND. IACT(K) .LE. M) THEN
                  if (temp > 0.0D0 .and. Iact(k) <= M) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      tempa = Vmultc(k) / temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              IF (RATIO .LT. 0.0 .OR. TEMPA .LT. RATIO) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-15: IOUT is never used
!                  IOUT=K
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      if (ratio < 0.0D0 .or. tempa < ratio) ratio = temp&
     &a
                  end if
                  if (k >= 2) then
                      kw = Iact(k)
                      do i = 1, N
                          Dxnew(i) = Dxnew(i) - temp * A(i, kw)
                      end do
                  end if
                  Vmultd(k) = temp
              else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          VMULTD(K)=0.0
                  Vmultd(k) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              end if
              k = k - 1
              if (k <= 0) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (RATIO .LT. 0.0) GOTO 490
                  if (ratio < 0.0D0) goto 700
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Revise the Lagrange multipliers and reorder the active constraints so
!     that the one to be replaced is at the end of the list. Also calculate the
!     new value of ZDOTA(NACT) and branch if it is not acceptable.
!
                  do k = 1, nact
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  160 VMULTC(K)=AMAX1(0.0,VMULTC(K)-RATIO*VMULTD(K))
                      Vmultc(k) = DMAX1(0.0D0, Vmultc(k) - ratio * Vmult&
     &d(k))
                  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  if (icon < nact) then
                      isave = Iact(icon)
                      vsave = Vmultc(icon)
                      k = icon
                      do
                          kp = k + 1
                          kw = Iact(kp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
                          sp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          do i = 1, N
                              sp = sp + Z(i, k) * A(i, kw)
                          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
                          temp = DSQRT(sp * sp + Zdota(kp)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          alpha = Zdota(kp) / temp
                          beta = sp / temp
                          Zdota(kp) = alpha * Zdota(k)
                          Zdota(k) = temp
                          do i = 1, N
                              temp = alpha * Z(i, kp) + beta * Z(i, k)
                              Z(i, kp) = alpha * Z(i, k) - beta * Z(i, k&
     &p)
                              Z(i, k) = temp
                          end do
                          Iact(k) = kw
                          Vmultc(k) = Vmultc(kp)
                          k = kp
                          if (k >= nact) then
                              Iact(k) = isave
                              Vmultc(k) = vsave
                              exit
                          end if
                      end do
                  end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
                  temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do i = 1, N
                      temp = temp + Z(i, nact) * A(i, kk)
                  end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TEMP .EQ. 0.0) GOTO 490
                  if (temp == 0.0D0) goto 700
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  Zdota(nact) = temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VMULTC(ICON)=0.0
                  Vmultc(icon) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  Vmultc(nact) = ratio
                  exit
              end if
          end do
      end if
!
!     Update IACT and ensure that the objective function continues to be
!     treated as the last active constraint when MCON>M.
!
300   Iact(icon) = Iact(nact)
      Iact(nact) = kk
      if (mcon > M .and. kk /= mcon) then
          k = nact - 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SP=0.0
          sp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i = 1, N
              sp = sp + Z(i, k) * A(i, kk)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=SQRT(SP*SP+ZDOTA(NACT)**2)
          temp = DSQRT(sp * sp + Zdota(nact)**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          alpha = Zdota(nact) / temp
          beta = sp / temp
          Zdota(nact) = alpha * Zdota(k)
          Zdota(k) = temp
          do i = 1, N
              temp = alpha * Z(i, nact) + beta * Z(i, k)
              Z(i, nact) = alpha * Z(i, k) - beta * Z(i, nact)
              Z(i, k) = temp
          end do
          Iact(nact) = Iact(k)
          Iact(k) = kk
          temp = Vmultc(k)
          Vmultc(k) = Vmultc(nact)
          Vmultc(nact) = temp
      end if
!
!     If stage one is in progress, then set SDIRN to the direction of the next
!     change to the current vector of variables.
!
      if (mcon <= M) then
          kk = Iact(nact)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
          temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i = 1, N
              temp = temp + Sdirn(i) * A(i, kk)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=TEMP-1.0
          temp = temp - 1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          temp = temp / Zdota(nact)
          do i = 1, N
              Sdirn(i) = Sdirn(i) - temp * Z(i, nact)
          end do
          goto 500
      end if
!
!     Pick the next search direction of stage two.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  320 TEMP=1.0/ZDOTA(NACT)
400   temp = 1.0D0 / Zdota(nact)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, N
          Sdirn(i) = temp * Z(i, nact)
      end do
!
!     Calculate the step to the boundary of the trust region or take the step
!     that reduces RESMAX to zero. The two statements below that include the
!     factor 1.0E-6 prevent some harmless underflows that occurred in a test
!     calculation. Further, we skip the step if it could be zero within a
!     reasonable tolerance for computer rounding errors.
!
500   dd = Rho * Rho
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SD=0.0
!      SS=0.0
      sd = 0.0D0
      ss = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (ABS(DX(I)) .GE. 1.0E-6*RHO) DD=DD-DX(I)**2
          if (DABS(Dx(i)) >= 1.0D-6 * Rho) dd = dd - Dx(i)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          sd = sd + Dx(i) * Sdirn(i)
          ss = ss + Sdirn(i)**2
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (DD .LE. 0.0) GOTO 490
!      TEMP=SQRT(SS*DD)
!      IF (ABS(SD) .GE. 1.0E-6*TEMP) TEMP=SQRT(SS*DD+SD*SD)
      if (dd <= 0.0D0) goto 700
      temp = DSQRT(ss * dd)
      if (DABS(sd) >= 1.0D-6 * temp) temp = DSQRT(ss * dd + sd * sd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      stpful = dd / (temp + sd)
      step = stpful
      if (mcon == M) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          ACCA=STEP+0.1*RESMAX
!          ACCB=STEP+0.2*RESMAX
          acca = step + 0.1D0 * resmax
          accb = step + 0.2D0 * resmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (step >= acca .or. acca >= accb) goto 600
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          STEP=AMIN1(STEP,RESMAX)
          step = DMIN1(step, resmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
!
!     Set DXNEW to the new variables if STEP is the steplength, and reduce
!     RESMAX to the corresponding maximum residual if stage one is being done.
!     Because DXNEW will be changed during the calculation of some Lagrange
!     multipliers, it will be restored to the following value later.
!
      do i = 1, N
          Dxnew(i) = Dx(i) + step * Sdirn(i)
      end do
      if (mcon == M) then
          resold = resmax
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RESMAX=0.0
          resmax = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k = 1, nact
              kk = Iact(k)
              temp = B(kk)
              do i = 1, N
                  temp = temp - A(i, kk) * Dxnew(i)
              end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RESMAX=AMAX1(RESMAX,TEMP)
              resmax = DMAX1(resmax, temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          end do
      end if
!
!     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
!     device is included to force VMULTD(K)=0.0 if deviations from this value
!     can be attributed to computer rounding errors. First calculate the new
!     Lagrange multipliers.
!
      k = nact
      do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  390 ZDOTW=0.0
!      ZDWABS=0.0
          zdotw = 0.0D0
          zdwabs = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i = 1, N
              temp = Z(i, k) * Dxnew(i)
              zdotw = zdotw + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  400 ZDWABS=ZDWABS+ABS(TEMP)
!      ACCA=ZDWABS+0.1*ABS(ZDOTW)
!      ACCB=ZDWABS+0.2*ABS(ZDOTW)
!      IF (ZDWABS .GE. ACCA .OR. ACCA .GE. ACCB) ZDOTW=0.0
              zdwabs = zdwabs + DABS(temp)
          end do
          acca = zdwabs + 0.1D0 * DABS(zdotw)
          accb = zdwabs + 0.2D0 * DABS(zdotw)
          if (zdwabs >= acca .or. acca >= accb) zdotw = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          Vmultd(k) = zdotw / Zdota(k)
          if (k >= 2) then
              kk = Iact(k)
              do i = 1, N
                  Dxnew(i) = Dxnew(i) - Vmultd(k) * A(i, kk)
              end do
              k = k - 1
              cycle
          end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (MCON .GT. M) VMULTD(NACT)=AMAX1(0.0,VMULTD(NACT))
          if (mcon > M) Vmultd(nact) = DMAX1(0.0D0, Vmultd(nact))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Complete VMULTC by finding the new constraint residuals.
!
          do i = 1, N
              Dxnew(i) = Dx(i) + step * Sdirn(i)
          end do
          if (mcon > nact) then
              kl = nact + 1
              do k = kl, mcon
                  kk = Iact(k)
                  sum = resmax - B(kk)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SUMABS=RESMAX+ABS(B(KK))
                  sumabs = resmax + DABS(B(kk))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do i = 1, N
                      temp = A(i, kk) * Dxnew(i)
                      sum = sum + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  430     SUMABS=SUMABS+ABS(TEMP)
!          ACCA=SUMABS+0.1*ABS(SUM)
!          ACCB=SUMABS+0.2*ABS(SUM)
!          IF (SUMABS .GE. ACCA .OR. ACCA .GE. ACCB) SUM=0.0
                      sumabs = sumabs + DABS(temp)
                  end do
                  acca = sumabs + 0.1D0 * DABS(sum)
                  accb = sumabs + 0.2D0 * DABS(sum)
                  if (sumabs >= acca .or. acca >= accb) sum = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  Vmultd(k) = sum
              end do
          end if
!
!     Calculate the fraction of the step from DX to DXNEW that will be taken.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=1.0
          ratio = 1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          icon = 0
          do k = 1, mcon
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (VMULTD(K) .LT. 0.0) THEN
              if (Vmultd(k) < 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  temp = Vmultc(k) / (Vmultc(k) - Vmultd(k))
                  if (temp < ratio) then
                      ratio = temp
                      icon = k
                  end if
              end if
          end do
!
!     Update DX, VMULTC and RESMAX.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=1.0-RATIO
          temp = 1.0D0 - ratio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i = 1, N
              Dx(i) = temp * Dx(i) + ratio * Dxnew(i)
          end do
          do k = 1, mcon
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  470 VMULTC(K)=AMAX1(0.0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
              Vmultc(k) = DMAX1(0.0D0, temp * Vmultc(k) + ratio * Vmultd&
     &(k))
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (mcon == M) resmax = resold + ratio * (resmax - resold)
!
!     If the full step is not acceptable then begin another iteration.
!     Otherwise switch to stage two or end the calculation.
!
          if (icon > 0) goto 200
          if (step /= stpful) exit
          goto 99999
      end do
600   mcon = M + 1
      icon = mcon
      Iact(mcon) = mcon
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VMULTC(MCON)=0.0
      Vmultc(mcon) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      goto 100
!
!     We employ any freedom that may be available to reduce the objective
!     function before returning a DX whose length is less than RHO.
!
700   if (mcon == M) goto 600
      Ifull = 0
99999 end subroutine TRSTLP