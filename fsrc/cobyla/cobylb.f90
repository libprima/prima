!*==cobylb.f90  processed by SPAG 7.50RE at 17:02 on 30 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SUBROUTINE COBYLB (N,M,MPP,X,RHOBEG,RHOEND,IPRINT,MAXFUN,
!     1  CON,SIM,SIMI,DATMAT,A,VSIG,VETA,SIGBAR,DX,W,IACT)
      subroutine COBYLB(N, M, Mpp, X, Rhobeg, Rhoend, Iprint, Maxfun, Con, Sim,  &
     &                  Simi, Datmat, A, Vsig, Veta, Sigbar, Dx, W, Iact, F, Info,&
     &                  Ftarget, Resmax)
      implicit none
!*--COBYLB14
!*++
!*++ PARAMETER definitions rewritten by SPAG
!*++
      integer, parameter :: NSMAX = 1000
      real*8, parameter :: CTOL = epsilon(1.0D0)
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      integer :: N
      integer :: M
      integer, intent(IN) :: Mpp
      real*8, intent(INOUT), dimension(*) :: X
      real*8, intent(IN) :: Rhobeg
      real*8, intent(IN) :: Rhoend
      integer, intent(IN) :: Iprint
      integer, intent(INOUT) :: Maxfun
      real*8, intent(INOUT), dimension(*) :: Con
      real*8, intent(INOUT), dimension(N, *) :: Sim
      real*8, intent(INOUT), dimension(N, *) :: Simi
      real*8, intent(INOUT), dimension(Mpp, *) :: Datmat
      real*8, intent(INOUT), dimension(N, *) :: A
      real*8, intent(INOUT), dimension(*) :: Vsig
      real*8, intent(INOUT), dimension(*) :: Veta
      real*8, intent(INOUT), dimension(*) :: Sigbar
      real*8, intent(INOUT), dimension(*) :: Dx
      real*8, intent(INOUT), dimension(*) :: W
      integer, dimension(*) :: Iact
      real*8, intent(INOUT) :: F
      integer, intent(INOUT) :: Info
      real*8, intent(IN) :: Ftarget
      real*8, intent(INOUT) :: Resmax
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      real*8 :: almost_infinity, alpha, barmu, beta, cmax, cmin,    &
     &        csum, cvmaxm, cvmaxp, delta, denom, dxsign, edgmax,&
     &        error, gamma, hugenum, pareta, parmu, parsig, phi, &
     &        phimin, prerec, prerem, ratio, resnew, resref, rho,&
     &        sum, temp, tempa, trured, vmnew, vmold, weta, wsig
      logical :: better
      real*8, dimension(Mpp) :: consav, datdrop
      real*8, dimension(Mpp, NSMAX) :: datsav
      integer :: i, ibrnch, idxnew, iflag, ifull, iptem, iptemp, &
     &           isdirn, ivmc, ivmd, iz, izdota, j, jdrop, k,   &
     &           l, mp, nbest, nfvals, np, nsav
      real*8, dimension(N) :: xdrop
      real*8, dimension(N, NSMAX) :: xsav
!*++
!*++ End of declarations rewritten by SPAG
!*++
! NSMAX is the maximal number of "dropped X" to save (see comments below
! line number 480)
! CTOL is the tolerance for consraint violation. A point X is considered
! to be feasible if its constraint violation (RESMAX) is less than CTOL.
! EPSILON(1.0D0) returns the machine epsilon corresponding to 1.0D0,
! which is expected to be about 2.0D-16.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  A(N,*),VSIG(*),VETA(*),SIGBAR(*),DX(*),W(*),IACT(*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set the initial values of some parameters. The last column of SIM holds
!     the optimal vertex of the current simplex, and the preceding N columns
!     hold the displacements from the optimal vertex to the other vertices.
!     Further, SIMI holds the inverse of the matrix that is contained in the
!     first N columns of SIM.
!
      Info = 2147483647
      iptem = MIN0(N, 5)
      iptemp = iptem + 1
      np = N + 1
      mp = M + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ALPHA=0.25
!      BETA=2.1
!      GAMMA=0.5
!      DELTA=1.1
      alpha = 0.25D0
      beta = 2.1D0
      gamma = 0.5D0
      delta = 1.1D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rho = Rhobeg
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      PARMU=0.0
      parmu = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (Iprint >= 2) print 99001, rho
99001 format(/3X, 'The initial value of RHO is', 1PE13.6, 2X,             &
     &        'and PARMU is set to zero.')
      nfvals = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=1.0/RHO
      temp = 1.0D0 / rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, N
          Sim(i, np) = X(i)
          do j = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SIM(I,J)=0.0
!   20 SIMI(I,J)=0.0
              Sim(i, j) = 0.0D0
              Simi(i, j) = 0.0D0
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          Sim(i, i) = rho
          Simi(i, i) = temp
      end do
      jdrop = np
      ibrnch = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      nsav = 0
      hugenum = huge(0.0D0)
      datsav = hugenum
      almost_infinity = huge(0.0D0) / 2.0D0
100   do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Make the next call of the user-supplied subroutine CALCFC. These
!     instructions are also used for calling CALCFC during the iterations of
!     the algorithm.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   40 IF (NFVALS .GE. MAXFUN .AND. NFVALS .GT. 0) THEN
!          IF (IPRINT .GE. 1) PRINT 50
!   50     FORMAT (/3X,'Return from subroutine COBYLA because the ',
!     1      'MAXFUN limit has been reached.')
!          GOTO 600
!      END IF
!      NFVALS=NFVALS+1
!      CALL CALCFC (N,M,X,F,CON)
!
!     By Zaikun (02-06-2019):
          do i = 1, N
              if (X(i) /= X(i)) then
                  F = X(i)
                  ! Set F to NaN
                  Info = -1
                  goto 500
              end if
          end do

          call CALCFC(N, M, X, F, Con)
          nfvals = nfvals + 1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RESMAX=0.0
          Resmax = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (M > 0) then
              do k = 1, M
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   60     RESMAX=AMAX1(RESMAX,-CON(K))
                  Resmax = DMAX1(Resmax, -Con(k))
              end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          end if


          if (nfvals == Iprint - 1 .or. Iprint == 3) then
              print 99009, nfvals, F, Resmax, (X(i), i=1, iptem)
              if (iptem < N) print 99010, (X(i), i=iptemp, N)
          end if
          Con(mp) = F
          Con(Mpp) = Resmax

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! By Zaikun 20190819:
! CONSAV always containts the containt value of the current x.
! CON, however, will be changed during the calculation (see the lines
! above line number 220).
          do k = 1, Mpp
              consav(k) = Con(k)
          end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom/Zaikun (on 04-06-2019/10-06-2019):
!     CSUM containts the sum of the absolute value of the constraints to
!     check whether it contains a NaN value.
          csum = 0.0D0
          do k = 1, M
              csum = csum + DABS(Con(k))
          end do
          if (csum /= csum) then
              Resmax = csum
              ! Set RESMAX to NaN
              Info = -2
              goto 500
          end if
!     If the objective function value or the constraints contain a NaN or an
!     infinite value, the algorithm stops.
          if (F /= F .or. F > almost_infinity) then
              Info = -2
              goto 500
          end if
!     If the objective function achieves the target value at a feasible
!     point, then exit.
!      IF (F .LE. FTARGET .AND. RESMAX .LE. 0.0D0) THEN
          if (F <= Ftarget .and. Resmax < CTOL) then
!         The feasibility is guarantee because RESMAX .LE. CTOL
              Info = 1
              goto 600
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Zaikun (on 06-06-2019)
!     The following code was placed before "CALL CALCFC (N,M,X,F,CON)".
!     This led to a bug, because F may not equal F(X) if the subroutine
!     exits due to NFVALS .GE. MAXFUN (X is updated but F is not evaluated
!     at X). Similar thing can be said about RESMAX.
          if (nfvals >= Maxfun .and. nfvals > 0) then
              if (Iprint >= 1) print 99002
99002         format(/3X, 'Return from subroutine COBYLA because the ',   &
       &              'MAXFUN limit has been reached.')
              Info = 3
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IF (IBRNCH == 1) GOTO 440
          if (ibrnch == 1 .and. Info /= 3) then
              vmold = Datmat(mp, np) + parmu * Datmat(Mpp, np)
              vmnew = F + parmu * Resmax
              trured = vmold - vmnew
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (PARMU .EQ. 0.0 .AND. F .EQ. DATMAT(MP,NP)) THEN
              if (parmu == 0.0D0 .and. F == Datmat(mp, np)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  prerem = prerec
                  trured = Datmat(Mpp, np) - Resmax
              end if
!
!     Begin the operations that decide whether x(*) should replace one of the
!     vertices of the current simplex, the change being mandatory if TRURED is
!     positive. Firstly, JDROP is set to the index of the vertex that is to be
!     replaced.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=0.0
!      IF (TRURED .LE. 0.0) RATIO=1.0
              ratio = 0.0D0
              if (trured <= 0.0D0) ratio = 1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              jdrop = 0
              do j = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
                  temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do i = 1, N
                      temp = temp + Simi(j, i) * Dx(i)
                  end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=ABS(TEMP)
                  temp = DABS(temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  if (temp > ratio) then
                      jdrop = j
                      ratio = temp
                  end if
                  Sigbar(j) = temp * Vsig(j)
              end do
!
!     Calculate the value of ell.
!
              edgmax = delta * rho
              l = 0
              do j = 1, N
                  if (Sigbar(j) >= parsig .or. Sigbar(j) >= Vsig(j)) then
                      temp = Veta(j)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TRURED .GT. 0.0) THEN
!              TEMP=0.0
                      if (trured > 0.0D0) then
                          temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          do i = 1, N
                              temp = temp + (Dx(i) - Sim(i, j))**2
                          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              TEMP=SQRT(TEMP)
                          temp = DSQRT(temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      end if
                      if (temp > edgmax) then
                          l = j
                          edgmax = temp
                      end if
                  end if
              end do
              if (l > 0) jdrop = l
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820:
! When JDROP=0, the algorithm decides not to include the trust-region
! trial point X into the simplex, because X is not good enough according
! to the merit function PHI = F + PARMU*RESMAX. In this case, X will
! simply be discarded in the original code. However, this decision
! depends on the value of PARMU. When PARMU is updated later, the
! discarded X might turn out better, sometimes even better than
! SIM(:, NP), which is supposed to be the best point in the simplex.
! For this reason, we save the to-be-discarded X in XSAV and compare
! them with SIM(:, NP) right before exiting. If a vector in XSAV turns
! out better than SIM(:, NP), we replace SIM(:, NP) by this vector
!
! When JDROP > 0, SIM(:, JDROP) will be removed from the simplex
! according to PHI with the current PARMU. Similar to X, SIM(:, JDROP)
! may turn out better when PARMU is updated. Therefore, XSAV also takes
! SIM(:, JDROP) into account.
!
! We save at most NSMAX to-be-discarded X.
!
              if (jdrop == 0) then
                  do i = 1, N
                      xdrop(i) = X(i)
                  end do
                  do k = 1, Mpp
                      datdrop(k) = consav(k)
                  end do
              else
                  ! JDROP < NP is guaranteed
                  do i = 1, N
                      xdrop(i) = Sim(i, np) + Sim(i, jdrop)
                  end do
                  do k = 1, Mpp
                      datdrop(k) = Datmat(k, jdrop)
                  end do
              end if
              call SAVEX(xdrop(1:N), datdrop(1:Mpp), xsav(1:N, 1:NSMAX),     &
       &                 datsav(1:Mpp, 1:NSMAX), N, M, nsav, NSMAX, CTOL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (jdrop == 0) goto 400
!
!     Revise the simplex by updating the elements of SIM, SIMI and DATMAT.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
              temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do i = 1, N
                  Sim(i, jdrop) = Dx(i)
                  temp = temp + Simi(jdrop, i) * Dx(i)
              end do
              do i = 1, N
                  Simi(jdrop, i) = Simi(jdrop, i) / temp
              end do
              do j = 1, N
                  if (j /= jdrop) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=0.0
                      temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      do i = 1, N
                          temp = temp + Simi(j, i) * Dx(i)
                      end do
                      do i = 1, N
                          Simi(j, i) = Simi(j, i) - temp * Simi(jdrop, i)
                      end do
                  end if
              end do
              do k = 1, Mpp
                  Datmat(k, jdrop) = Con(k)
              end do
!
!     Branch back for further iterations with the current RHO.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TRURED .GT. 0.0 .AND. TRURED .GE. 0.1*PREREM) GOTO 140
              if (trured <= 0.0D0 .or. trured < 0.1D0 * prerem) goto 400
              exit
          else
!
!     Set the recently calculated function values in a column of DATMAT. This
!     array has a column for each vertex of the current simplex, the entries of
!     each column being the values of the constraint functions (if any)
!     followed by the objective function and the greatest constraint violation
!     at the vertex.
!
              do k = 1, Mpp
                  Datmat(k, jdrop) = Con(k)
              end do

              if (nfvals <= np) then
                  ! IF we do not go to 130 but continue to below, then NFVALS <= NP.
                  ! Thus NFVALS may be NP = N+1 > N.
!
!     Exchange the new vertex of the initial simplex with the optimal vertex if
!     necessary. Then, if the initial simplex is not complete, pick its next
!     vertex and calculate the function values there.
!
                  if (jdrop <= N) then
                      if (Datmat(mp, np) <= F) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! When NFVALS<=N, it is necessary to update X(JDROP) because next X will
! be calculated based on the current one (see the code below line number
! 120). The purpose of this update is to make this X hold the variable that
! has the smallest function value up to now. The next X is defined as
! a perturbation of this optimal X, which is reasonable.
! However, this update leads to inconsistency between X and [F, RESMAX],
! meaning that X is not necessarily the X corresponding to [F, RESMAX].
! This can cause COBYLA return inconsistent [X, F, RESMAX].
! Fortunately, this is not a problem if NFVALS <= N. Because, if COBYLA
! returns with a NFVALS <= N, then X contains NaN or F = NaN or nearly
! Inf or the constraint contains NaN, all of which would lead to an
! immediate jump to line 600 without coming here. Therefore, if we
! restrict the update to only the case with NFVALS <= N, ther will be no
! inconsistency at return.
! With the original code, if COBYLA returns with NFVALS = NP (this can
! happen if the initial TR problem constantly produces too short steps
! so that RHO is reduced to RHOEND without producing any acceptable trial
! step; see the code below line number 380), then, when the code arrives
! at line number 600, X and [F, RESMAX] may be inconsistent. However,
! recall that the inconsistency occurs only if SIM(:, NP) has a lower
! function value than X (before updated). Therefore, as long as the code
! takes SIM(:, NP) instead of X, no harm would be done. It is the case
! likely the in the original code, because if COBYLA returns with
! NFVALS=NP, then PARMU=0, and hence SIM(:, NP) will be returned.
!
!              X(JDROP)=SIM(JDROP,NP)
                          if (nfvals <= N) X(jdrop) = Sim(jdrop, np)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      else
                          Sim(jdrop, np) = X(jdrop)
                          do k = 1, Mpp
                              Datmat(k, jdrop) = Datmat(k, np)
                              Datmat(k, np) = Con(k)
                          end do
                          do k = 1, jdrop
                              Sim(jdrop, k) = -rho
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              TEMP=0.0
                              temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              do i = k, jdrop
                                  temp = temp - Simi(i, k)
                              end do
                              Simi(jdrop, k) = temp
                          end do
                      end if
                  end if

! 120
                  if (nfvals <= N) then
                      jdrop = nfvals
                      X(jdrop) = X(jdrop) + rho
                      cycle
                  end if
              end if
              ibrnch = 1
              exit
          end if
      end do
200   do
!
!     Identify the optimal vertex of the current simplex.
!
          phimin = Datmat(mp, np) + parmu * Datmat(Mpp, np)
          nbest = np
          do j = 1, N
              temp = Datmat(mp, j) + parmu * Datmat(Mpp, j)
              if (temp < phimin) then
                  nbest = j
                  phimin = temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ELSE IF (TEMP .EQ. PHIMIN .AND. PARMU .EQ. 0.0) THEN
              elseif (temp == phimin .and. parmu == 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  if (Datmat(Mpp, j) < Datmat(Mpp, nbest)) nbest = j
              end if
          end do
!
!     Switch the best vertex into pole position if it is not there already,
!     and also update SIM, SIMI and DATMAT.
!
          if (nbest <= N) then
              do i = 1, Mpp
                  temp = Datmat(i, np)
                  Datmat(i, np) = Datmat(i, nbest)
                  Datmat(i, nbest) = temp
              end do
              do i = 1, N
                  temp = Sim(i, nbest)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SIM(I,NBEST)=0.0
                  Sim(i, nbest) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  Sim(i, np) = Sim(i, np) + temp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMPA=0.0
                  tempa = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do k = 1, N
                      Sim(i, k) = Sim(i, k) - temp
                      tempa = tempa - Simi(k, i)
                  end do
                  Simi(nbest, i) = tempa
              end do
          end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          if (Info == 3) goto 500
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Make an error return if SIGI is a poor approximation to the inverse of
!     the leading N by N submatrix of SIG.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ERROR=0.0
          error = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i = 1, N
              do j = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
!      IF (I .EQ. J) TEMP=TEMP-1.0
                  temp = 0.0D0
                  if (i == j) temp = temp - 1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do k = 1, N
                      temp = temp + Simi(i, k) * Sim(k, j)
                  end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  200 ERROR=AMAX1(ERROR,ABS(TEMP))
!      IF (ERROR .GT. 0.1) THEN
                  error = DMAX1(error, DABS(temp))
              end do
          end do
          if (error > 0.1D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (Iprint >= 1) print 99003
99003         format(/3X, 'Return from subroutine COBYLA because ',       &
       &              'rounding errors are becoming damaging.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
              Info = 7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              goto 500
          end if
!
!     Calculate the coefficients of the linear approximations to the objective
!     and constraint functions, placing minus the objective function gradient
!     after the constraint gradients in the array A. The vector W is used for
!     working space.
!
! 220
          do k = 1, mp
              Con(k) = -Datmat(k, np)
              do j = 1, N
                  W(j) = Datmat(k, j) + Con(k)
              end do
              do i = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
                  temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do j = 1, N
                      temp = temp + W(j) * Simi(j, i)
                  end do
                  if (k == mp) temp = -temp
                  A(i, k) = temp
              end do
          end do
!
!     Calculate the values of sigma and eta, and set IFLAG=0 if the current
!     simplex is not acceptable.
!
          iflag = 1
          parsig = alpha * rho
          pareta = beta * rho
          do j = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      WSIG=0.0
!      WETA=0.0
              wsig = 0.0D0
              weta = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do i = 1, N
                  wsig = wsig + Simi(j, i)**2
                  weta = weta + Sim(i, j)**2
              end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VSIG(J)=1.0/SQRT(WSIG)
!      VETA(J)=SQRT(WETA)
              Vsig(j) = 1.0D0 / DSQRT(wsig)
              Veta(j) = DSQRT(weta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (Vsig(j) < parsig .or. Veta(j) > pareta) iflag = 0
          end do
!
!     If a new vertex is needed to improve acceptability, then decide which
!     vertex to drop from the simplex.
!
          if (ibrnch == 1 .or. iflag == 1) then
!
!     Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO.
!
              iz = 1
              izdota = iz + N * N
              ivmc = izdota + N
              isdirn = ivmc + mp
              idxnew = isdirn + N
              ivmd = idxnew + N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRSTLP is not predictable, and Segmentation Fault or infinite
! cycling may happen. This is because any equality/inequality comparison
! involving NaN returns FALSE, which can lead to unintended behavior of
! the code, including uninitialized indices.
              do j = 1, N
                  do i = 1, N
                      if (Simi(i, j) /= Simi(i, j)) then
                          if (Iprint >= 1) print 99004
99004                     format(/3X,                                       &
          &                       'Return from subroutine COBYLA because ',  &
          &                       'rounding errors are becoming damaging.')
                          Info = 7
                          goto 500
                      end if
                  end do
              end do
              do j = 1, mp
                  do i = 1, N
                      if (A(i, j) /= A(i, j)) then
                          Info = -3
                          goto 500
                      end if
                  end do
              end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              call TRSTLP(N, M, A, Con, rho, Dx, ifull, Iact, W(iz), W(izdota),    &
       &                  W(ivmc), W(isdirn), W(idxnew), W(ivmd))
              if (ifull == 0) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=0.0
                  temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do i = 1, N
                      temp = temp + Dx(i)**2
                  end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .LT. 0.25*RHO*RHO) THEN
                  if (temp < 0.25D0 * rho * rho) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      ibrnch = 1
                      exit
                  end if
              end if
!
!     Predict the change to F and the new maximum constraint violation if the
!     variables are altered from x(0) to x(0)+DX.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RESNEW=0.0
!      CON(MP)=0.0
              resnew = 0.0D0
              Con(mp) = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do k = 1, mp
                  sum = Con(k)
                  do i = 1, N
                      sum = sum - A(i, k) * Dx(i)
                  end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (K .LT. MP) RESNEW=AMAX1(RESNEW,SUM)
                  if (k < mp) resnew = DMAX1(resnew, sum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              end do
!
!     Increase PARMU if necessary and branch back if this change alters the
!     optimal vertex. Otherwise PREREM and PREREC will be set to the predicted
!     reductions in the merit function and the maximum constraint violation
!     respectively.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      BARMU=0.0
              barmu = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              prerec = Datmat(Mpp, np) - resnew
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (PREREC .GT. 0.0) BARMU=SUM/PREREC
!      IF (PARMU .LT. 1.5*BARMU) THEN
!          PARMU=2.0*BARMU
              if (prerec > 0.0D0) barmu = sum / prerec
              if (parmu < 1.5D0 * barmu) then
                  parmu = 2.0D0 * barmu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  if (Iprint >= 2) print 99005, parmu
99005             format(/3X, 'Increase in PARMU to', 1PE13.6)
                  phi = Datmat(mp, np) + parmu * Datmat(Mpp, np)
                  do j = 1, N
                      temp = Datmat(mp, j) + parmu * Datmat(Mpp, j)
                      if (temp < phi) goto 300
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .EQ. PHI .AND. PARMU .EQ. 0.0) THEN
                      if (temp == phi .and. parmu == 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          if (Datmat(Mpp, j) < Datmat(Mpp, np)) goto 300
                      end if
                  end do
              end if
              prerem = parmu * prerec - sum
!
!     Calculate the constraint and objective functions at x(*). Then find the
!     actual reduction in the merit function.
!
              do i = 1, N
                  X(i) = Sim(i, np) + Dx(i)
              end do
              ibrnch = 1
          else
              jdrop = 0
              temp = pareta
              do j = 1, N
                  if (Veta(j) > temp) then
                      jdrop = j
                      temp = Veta(j)
                  end if
              end do
              if (jdrop == 0) then
                  do j = 1, N
                      if (Vsig(j) < temp) then
                          jdrop = j
                          temp = Vsig(j)
                      end if
                  end do
              end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190822: If VETA or VSIG become NaN due to rounding errors,
! JDROP may end up being 0. If we continue, then a Segmentation Fault
! will happen because we will read SIM(:, JDROP) and VSIG(JDROP).
              if (jdrop == 0) then
                  if (Iprint >= 1) print 99006
99006             format(/3X, 'Return from subroutine COBYLA because ',    &
        &                 'rounding errors are becoming damaging.')
                  Info = 7
                  goto 500
              end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820: See the comments below line number 480
              do i = 1, N
                  xdrop(i) = Sim(i, np) + Sim(i, jdrop)
                  ! JDROP<NP is guaranteed
              end do
              do k = 1, Mpp
                  datdrop(k) = Datmat(k, jdrop)
              end do
              call SAVEX(xdrop(1:N), datdrop(1:Mpp), xsav(1:N, 1:NSMAX),     &
       &                 datsav(1:Mpp, 1:NSMAX), N, M, nsav, NSMAX, CTOL)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the step to the new vertex and its sign.
!
              temp = gamma * rho * Vsig(jdrop)
              do i = 1, N
                  Dx(i) = temp * Simi(jdrop, i)
              end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      CVMAXP=0.0
!      CVMAXM=0.0
              cvmaxp = 0.0D0
              cvmaxm = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do k = 1, mp
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SUM=0.0
                  sum = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do i = 1, N
                      sum = sum + A(i, k) * Dx(i)
                  end do
                  if (k < mp) then
                      temp = Datmat(k, np)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          CVMAXP=AMAX1(CVMAXP,-SUM-TEMP)
!          CVMAXM=AMAX1(CVMAXM,SUM-TEMP)
                      cvmaxp = DMAX1(cvmaxp, -sum - temp)
                      cvmaxm = DMAX1(cvmaxm, sum - temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  end if
              end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      DXSIGN=1.0
!      IF (PARMU*(CVMAXP-CVMAXM) .GT. SUM+SUM) DXSIGN=-1.0
              dxsign = 1.0D0
              if (parmu * (cvmaxp - cvmaxm) > sum + sum) dxsign = -1.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Update the elements of SIM and SIMI, and set the next X.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
              temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do i = 1, N
                  Dx(i) = dxsign * Dx(i)
                  Sim(i, jdrop) = Dx(i)
                  temp = temp + Simi(jdrop, i) * Dx(i)
              end do
              do i = 1, N
                  Simi(jdrop, i) = Simi(jdrop, i) / temp
              end do
              do j = 1, N
                  if (j /= jdrop) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=0.0
                      temp = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      do i = 1, N
                          temp = temp + Simi(j, i) * Dx(i)
                      end do
                      do i = 1, N
                          Simi(j, i) = Simi(j, i) - temp * Simi(jdrop, i)
                      end do
                  end if
                  X(j) = Sim(j, np) + Dx(j)
              end do
          end if
          goto 100
300   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
400   if (iflag == 0) then
          ibrnch = 0
          goto 200
      end if
!
!     Otherwise reduce RHO if it is not at its least value and reset PARMU.
!
      if (rho > Rhoend) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RHO=0.5*RHO
!          IF (RHO .LE. 1.5*RHOEND) RHO=RHOEND
!          IF (PARMU .GT. 0.0) THEN
!              DENOM=0.0
          rho = 0.5D0 * rho
          if (rho <= 1.5D0 * Rhoend) rho = Rhoend
          if (parmu > 0.0D0) then
              denom = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do k = 1, mp
                  cmin = Datmat(k, np)
                  cmax = cmin
                  do i = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              CMIN=AMIN1(CMIN,DATMAT(K,I))
!  560         CMAX=AMAX1(CMAX,DATMAT(K,I))
!              IF (K .LE. M .AND. CMIN .LT. 0.5*CMAX) THEN
!                  TEMP=AMAX1(CMAX,0.0)-CMIN
!                  IF (DENOM .LE. 0.0) THEN
                      cmin = DMIN1(cmin, Datmat(k, i))
                      cmax = DMAX1(cmax, Datmat(k, i))
                  end do
                  if (k <= M .and. cmin < 0.5D0 * cmax) then
                      temp = DMAX1(cmax, 0.0D0) - cmin
                      if (denom <= 0.0D0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          denom = temp
                      else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                      DENOM=AMIN1(DENOM,TEMP)
                          denom = DMIN1(denom, temp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      end if
                  end if
              end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              IF (DENOM .EQ. 0.0) THEN
!                  PARMU=0.0
              if (denom == 0.0D0) then
                  parmu = 0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              elseif (cmax - cmin < parmu * denom) then
                  parmu = (cmax - cmin) / denom
              end if
          end if
          if (Iprint >= 2) print 99007, rho, parmu
99007     format(/3X, 'Reduction in RHO to', 1PE13.6, '  and PARMU =',     &
      &           1PE13.6)
          if (Iprint == 2) then
              print 99009, nfvals, Datmat(mp, np), Datmat(Mpp, np),     &
       &            (Sim(i, np), i=1, iptem)
              if (iptem < N) print 99010, (X(i), i=iptemp, N)
          end if
          goto 200
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      else
          Info = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
!
!     Return the best calculated values of the variables.
!
      if (Iprint >= 1) print 99008
99008 format(/3X, 'Normal return from subroutine COBYLA')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (IFULL .EQ. 1) GOTO 620
!  600 DO 610 I=1,N
!  610 X(I)=SIM(I,NP)
!      F=DATMAT(MP,NP)
!      RESMAX=DATMAT(MPP,NP)
!
!      Zaikun 01-06-2019:
!      Why go to 620 directly without setting X and F? This seems
!      INCORRECT, because it may lead to a return with X and F that
!      are not the best available.
!      The following code defines X as an "optimal" one among the
!      vectors:
!      DATSAV(:, 1:NSAV), (when NSAV >= 1),
!      SIM(:, NP), and SIM(:, 1:min(NP-1, NFVALS-2)) (when NFVALS>=2).
!      Here, X being "optimal" means
!      1. the constraint violation of X is at most RESREF
!      2. no other vector is better than X according to ISBETTER with
!      the current PARMU.
!
!      Note:
!      0. The last evaluated X and its function/constraint information
!      are saved in [X, CONSAV, F, RESMAX].
!      1. When NFVALS=1, SIM and DATMAT have not been initialized yet.
!      2. When 2<=NFVALS<=NP, the first evaluated X are saved in
!      SIM(:, NP), its function/constraint in DATMAT(:, NP), while the
!      other X are saved in SIM(:, NFVALS-1), its function/constraint
!      in DATMAT(:, NFVALS-1). However, when the code arrives at line 600,
!      [X, CON, F, RESMAX] may have not been saved into SIM(:, NFVALS-1)
!      and DATMAT(:, NFVALS-1) yet. That is why we check SIM up to
!      NFVALS-2 instead of NFVALS-1.
500   do k = 1, M
          Con(k) = consav(k)
      end do
      parmu = max(parmu, 1.0D2)
      if (nfvals >= 2) then ! See the comments above for why NFVALS>2
          call ISBETTER(F, Resmax, Datmat(mp, np), Datmat(Mpp, np), parmu, CTOL,&
      &                 better)
          if (better) then
              do i = 1, N
                  X(i) = Sim(i, np)
              end do
              F = Datmat(mp, np)
              Resmax = Datmat(Mpp, np)
              do k = 1, M
                  Con(k) = Datmat(k, np)
              end do
          end if
          resref = Resmax
          if (resref /= resref) resref = hugenum
          do j = 1, min(np - 1, nfvals - 2)
! See the comments above for why to check these J
              if (Datmat(Mpp, j) <= resref) then
                  call ISBETTER(F, Resmax, Datmat(mp, j), Datmat(Mpp, j), parmu, &
        &                       CTOL, better)
                  if (better) then
                      do i = 1, N
                          X(i) = Sim(i, j) + Sim(i, np)
                      end do
                      F = Datmat(mp, j)
                      Resmax = Datmat(Mpp, j)
                      do k = 1, M
                          Con(k) = Datmat(k, j)
                      end do
                  end if
              end if
          end do
      end if
      if (nsav >= 1) then ! Do the following only if NSAV >= 1.
!          DO J = 1, NSAV
          do j = nsav, 1, -1 ! We start with the most recent point
              if (datsav(Mpp, j) <= resref) then
                  call ISBETTER(F, Resmax, datsav(mp, j), datsav(Mpp, j), parmu, &
        &                       CTOL, better)
                  if (better) then
                      do i = 1, N
                          X(i) = xsav(i, j)
                      end do
                      F = datsav(mp, j)
                      Resmax = datsav(Mpp, j)
                      do k = 1, M
                          Con(k) = datsav(k, j)
                      end do
                  end if
              end if
          end do
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
600   if (Iprint >= 1) then
          print 99009, nfvals, F, Resmax, (X(i), i=1, iptem)
          if (iptem < N) print 99010, (X(i), i=iptemp, N)
      end if
      Maxfun = nfvals
99009 format(/3X, 'NFVALS =', I5, 3X, 'F =', 1PE13.6, 4X, 'MAXCV =',          &
     &        1PE13.6 / 3X, 'X =', 1PE13.6, 1P4E15.6)
99010 format(1PE19.6, 1P4E15.6)
      end subroutine COBYLB
!*==savex.f90  processed by SPAG 7.50RE at 17:02 on 30 May 2021

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820: See the comments below line number 480
      subroutine SAVEX(Xdrop, Datdrop, Xsav, Datsav, N, M, Nsav, Nsmax, Ctol)
! This subroutine saves XDROP in XSAV and DATDROP in DATSAV, unless
! XDROP is dominated by a vector in XSAV(:, 1:NSAV). If XDROP dominates
! some vectors in XSAV(:, 1:NSAV), then these vectors will be removed.
! If XDROP does not dominate any of XSAV(:, 1:NSAV) but NSAV=NSMAX,
! then we remove XSAV(:,1), which is the oldest vector in XSAV(:, 1:NSAV).
!
! When COBYLA calls this subroutine, XDROP is a vector to be "dropped",
! and  DATDROP contains its function/constraint inforation (in particular,
! DATDROP(MP) = F(XDROP), and DATDROP(MPP) = RESMAX(X)). XSAV and DATSAV
! save at most NSMAX vectors "dropped" by COBYLB and their function/constraint
! information. Only XSAV(:, 1:NSAV) and DATSAV(:, 1:NSAV) contains such
! vectors, while XSAV(:, NSAV+1:NSMAX) and DATSAV(:, NSAV+1:NSMAX) are
! not initialized yet.
!
! Note: X dominates Y if and only if the function/constraint of X is
! better than the function/constraint of Y accoring to the ISBETTER
! subroutine with PARMU = -1.0D0. Indeed, PARMU can be any negative
! number. This is because, due to the implementation of ISBETTER,
! X dominates Y (i.e., X is better than Y with PARMU < 0)
! ==> X is better than Y with any PARMU >= 0,
! ==> X is better than Y regardless of PARMU.
! Fot this reason, it is sufficient to save all the "dropped" X that are not
! dominated by any vectors in XSAV (as we do in this subroutine),
! because the other X are always worse than some vector in XSAV.
!
      implicit none
!*--SAVEX1042
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      real(kind(0.0D0)), intent(IN), dimension(N) :: Xdrop
      real(kind(0.0D0)), intent(IN), dimension(M + 2) :: Datdrop
      real(kind(0.0D0)), intent(INOUT), dimension(N, Nsmax) :: Xsav
      real(kind(0.0D0)), intent(INOUT), dimension(M + 2, Nsmax) :: Datsav
      integer, intent(IN) :: N
      integer, intent(IN) :: M
      integer, intent(INOUT) :: Nsav
      integer, intent(IN) :: Nsmax
      real(kind(0.0D0)), intent(IN) :: Ctol
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      logical :: better
      integer :: i, j, k, l, mp, mpp, nremove
      integer, dimension(Nsmax) :: iremove
      real(kind(0.0D0)) :: parmu
!*++
!*++ End of declarations rewritten by SPAG
!*++

      if (Nsmax <= 0) return ! Do nothing if NSMAX=0

      mp = M + 1
      mpp = M + 2
      parmu = -1.0D0 ! See the comments above for why PARMU = -1

! IREMOVE: indices of vectors to remove from XSAV
! NREMOVE: number of vectors to remove from XSAV
      do j = 1, Nsmax
! It is not enough to initialize IREMOVE(1:NSAV), because NSAV may be
! incremented by 1 latter, and then IREMOVE(NSAV+1) will be accessed.
          iremove(j) = -1
      end do
      nremove = 0
      do i = 1, Nsav
! If XDROP is dominated by XSAV(:, I), then return immediately,
! because XDROP should not be inluded into XSAV.
          call ISBETTER(Datdrop(mp), Datdrop(mpp), Datsav(mp, i),           &
      &                 Datsav(mpp, i), parmu, Ctol, better)
          if (better) return
! If XDROP dominates XSAV(:, I), then increment NREMOVE by 1 and save
! I as IREMOVE(NREMOVE).
          call ISBETTER(Datsav(mp, i), Datsav(mpp, i), Datdrop(mp),          &
      &                 Datdrop(mpp), parmu, Ctol, better)
          if (better) then
              nremove = nremove + 1
              iremove(nremove) = i
          end if
      end do

! The code did not return and NREMOVE=0 (no vector to remove from XSAV).
! If NSAV=NSMAX, then remove XSAV(:, 1); otherwise, increment NSAV by
! 1 and then "remove" XSAV(:, NSAV) (though there is no vector saved there)
      if (nremove == 0) then
          if (Nsav == Nsmax) then
              iremove(1) = 1
          else
              Nsav = Nsav + 1
              iremove(1) = Nsav
          end if
          nremove = 1
      end if

! Remove from XSAV the vectors indexed by IREMOVE
      j = 1 ! Index of IREMOVE
      k = 1 ! Index of the new XSAV
      do i = 1, Nsav
          ! Index of the old XSAV
          if (i == iremove(j)) then
              j = j + 1   ! Do nothing but incrementing J by 1
          else  ! Set XSAV(:, K) = XSAV(:, I)
              do l = 1, N
                  Xsav(l, k) = Xsav(l, i)
              end do
              do l = 1, mpp
                  Datsav(l, k) = Datsav(l, i)
              end do
              k = k + 1   ! Increment K by 1
          end if
      end do

! Set the number of vectors in the new XSAV
      Nsav = Nsav - nremove + 1

! Save XDROP in XSAV(:, NSAV) (with NSAV updated as above)
      if (Nsav >= 1 .and. Nsav <= Nsmax) then
          ! This inequlity is not guaranteed if NSMAX=0, where NSAV will
          ! be 0 and hence a Segmentation Fault when accessing
          ! XSAV(:, NSAV). Although we return immediately if NSMAX=0,
          ! we still check this inequlity to be safe.
          do l = 1, N
              Xsav(l, Nsav) = Xdrop(l)
          end do
          do l = 1, mpp
              Datsav(l, Nsav) = Datdrop(l)
          end do
      end if

      end subroutine SAVEX
!*==isbetter.f90  processed by SPAG 7.50RE at 17:02 on 30 May 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820:
      subroutine ISBETTER(F0, R0, F, R, Parmu, Ctol, Better)
! This subroutine compares whether (F, R) is (strictly) better than
! (F0, R0) in the sense of decreasing the merit function PHI = F + PARMU*R.
! It takes care of the cases where some of these values are NaN or Inf.
! At return, BETTER=true iff (F, R) is better than (F0, R0).
! Note:
! 1. A = Inf if and only if A .GT. HUGENUM (defined below);
! 2. A = NaN if and only if A .NE. A;
! 3. If A = NaN, then any comparison (except .NE.) with another number B
!    (can be Inf or NaN) returns false.
!
      implicit none
!*--ISBETTER1162
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
      real(kind(0.0D0)), intent(IN) :: F0
      real(kind(0.0D0)), intent(IN) :: R0
      real(kind(0.0D0)), intent(IN) :: F
      real(kind(0.0D0)), intent(IN) :: R
      real(kind(0.0D0)), intent(IN) :: Parmu
      real(kind(0.0D0)), intent(IN) :: Ctol
      logical, intent(OUT) :: Better
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
      logical :: f0infnan, finfnan, fle, flt, r0infnan, rinfnan,  &
     &           rle, rlt
      real(kind(0.0D0)) :: hugenum = huge(0.0D0)
!*++
!*++ End of declarations rewritten by SPAG
!*++

      Better = .false.

! As values of F0, R0, F, and R, we regard Inf and NaN being equivalent
! values (they are equally bad).
      f0infnan = (F0 /= F0) .or. (F0 > hugenum)     ! F0 = Inf or NaN?
      r0infnan = (R0 /= R0) .or. (R0 > hugenum)     ! R0 = Inf or NaN?
      finfnan = (F /= F) .or. (F > hugenum)     ! F = Inf or NaN?
      rinfnan = (R /= R) .or. (R > hugenum)     ! R  = Inf or NaN?

! When PARMU >= 0 and F + PARMU*R < F0 + PARMU*R0 and R < CTOL (feasible),
! then (F, R) is better than (F0, R0).
! Note that we should not set BETTER=FALSE even if this inequlity does not
! hold, because one or both of the two sides may be NaN.
      if (Parmu >= 0.0D0 .and. F + Parmu * R < F0 + Parmu * R0 .and. R < Ctol)      &
     &     Better = .true.

! If R < CTOL and F is not Inf or NaN while (R0 < CTOL) is false (may
! be because R0 is NaN), then (F, R) is better than (F0, R0). We prefer
! feasible points (i.e., constraint violation is less than CTOL) to
! insfeasible ones.
      if (R < Ctol .and. .not. (R0 < Ctol) .and. .not. finfnan)             &
     &     Better = .true.

! If F0 or R0 is Inf/NaN while neither F nor R is Inf/NaN, then (F, R)
! is better than (F0, R0).
      if ((f0infnan .or. r0infnan) .and. .not. (finfnan .or. rinfnan)) &
     &     Better = .true.

      flt = (f0infnan .and. (.not. finfnan)) .or. (F < F0)    ! F < F0?
      fle = (f0infnan .and. finfnan) .or. (F <= F0) .or. flt  ! F <= F0?
      rlt = (r0infnan .and. (.not. rinfnan)) .or. (R < R0)    ! R < R0?
      rle = (r0infnan .and. rinfnan) .or. (R <= R0) .or. rlt  ! R <= R0?

! If (F < F0 and R <= R0) or (F <= F0 and R < R0) in the sense defined
! above, the (F, R) is better than (F0, R0).
      if ((flt .and. rle) .or. (fle .and. rlt)) Better = .true.

! If one of F and R is -Inf and the other is not Inf/Nan, while neither
! F0 nor R0 is -Inf, then the (F, R) is better than (F0, R0).
      if ((F < -hugenum) .and. (.not. rinfnan) .and. (F0 >= -hugenum) .and. &
     &     (R0 >= -hugenum)) Better = .true.
      if ((R < -hugenum) .and. (.not. finfnan) .and. (F0 >= -hugenum) .and. &
     &     (R0 >= -hugenum)) Better = .true.
      end subroutine ISBETTER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
