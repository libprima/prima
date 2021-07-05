!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of cobylb.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 04-Jul-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      module cobylb_mod


      contains

      subroutine cobylb(m, x, rhobeg, rhoend, iprint, maxfun, con, f, in&
     &fo, ftarget, resmax)

! Generic modules
      use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, HUGENU&
     &M, DEBUGGING, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : calquad, inprod, matprod
      use logging_mod, only : logging

! Solver-specific modules
!use savex_mod, only : savex
!use isbetter_mod, only : isbetter
      use trustregion_mod, only : trstlp

      implicit none

! Inputs
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: m
      integer(IK), intent(in) :: maxfun
      real(RP), intent(in) :: ftarget
      real(RP), intent(in) :: rhobeg
      real(RP), intent(in) :: rhoend

! In-outputs
      real(RP), intent(out) :: con(:)
      ! m+2. Bad name; should be confr
      real(RP), intent(inout) :: x(:)
      ! n

! Outputs
      integer(IK), intent(out) :: info
      real(RP), intent(out) :: f


! Parameters
! NSMAX is the maximal number of "dropped X" to save (see comments below line number 480).
      integer(IK), parameter :: NSMAX = 1000
! CTOL is the tolerance for constraint violation. A point X is considered to be feasible if its
! constraint violation (RESMAX) is less than CTOL.
      real(RP), parameter :: CTOL = epsilon(1.0_RP)

! Local variables

      integer(IK) :: i
      integer(IK) :: iact(m + 1)
      integer(IK) :: ibrnch
      integer(IK) :: idxnew
      integer(IK) :: ifull
      integer(IK) :: iptem
      integer(IK) :: isdirn
      integer(IK) :: ivmc
      integer(IK) :: ivmd
      integer(IK) :: iz
      integer(IK) :: izdota
      integer(IK) :: j
      integer(IK) :: jdrop
      integer(IK) :: jmax
      integer(IK) :: jopt
      integer(IK) :: k
      integer(IK) :: l
      integer(IK) :: n
      integer(IK) :: nf
      integer(IK) :: nsav
      real(RP) :: A(size(x), m + 1)
      ! (n, )
! A(:, 1:m) contains the approximate gradient for the constraints, and A(:, m+1) is minus the
! approximate gradient for the objective function.
      real(RP) :: almost_infinity
      real(RP) :: alpha
      real(RP) :: barmu
      real(RP) :: beta
      real(RP) :: cmax
      real(RP) :: cmin
      real(RP) :: consav(m + 2)
      real(RP) :: cvmaxm
      real(RP) :: cvmaxp
      real(RP) :: datmat(m + 2, size(x) + 1)
      !(mpp, )
      real(RP) :: datsav(m + 2, nsmax)
      real(RP) :: delta
      real(RP) :: denom
      real(RP) :: dx(size(x))
      real(RP) :: dxsign
      real(RP) :: edgmax
      real(RP) :: error(size(x), size(x))
      real(RP) :: gamma
      real(RP) :: pareta
      real(RP) :: parmu
      real(RP) :: parsig
      real(RP) :: phi(size(x) + 1)
      real(RP) :: phitmp
      real(RP) :: phimin
      real(RP) :: prerec
      real(RP) :: prerem
      real(RP) :: ratio
      real(RP) :: resmax
      real(RP) :: resnew
      real(RP) :: resref
      real(RP) :: rho
      real(RP) :: sigbar(size(x))
      real(RP) :: sim(size(x), size(x) + 1)
      ! (n, )
      real(RP) :: simi(size(x), size(x))
      ! (n, )
      real(RP) :: summ
      real(RP) :: temp
      real(RP) :: simjopt(size(x))
      real(RP) :: trured
      real(RP) :: veta(size(x))
      real(RP) :: vmnew
      real(RP) :: vmold
      real(RP) :: vsig(size(x))
      real(RP) :: w(size(x) * (3 * size(x) + 2 * m + 11) + 4 * m + 6)
      real(RP) :: xsav(size(x), NSMAX)

      logical :: improve_geo

      character(len=SRNLEN), parameter :: srname = 'COBYLB'

      n = size(x)

! Set the initial values of some parameters. The last column of SIM holds the optimal vertex of the
! current simplex, and the preceding N columns hold the displacements from the optimal vertex to the
! other vertices.  Further, SIMI holds the inverse of the matrix that is contained in the first N
! columns of SIM.
      info = 2147483647
      iptem = min(n, 5)
      alpha = 0.25E0_RP
      beta = 2.1E0_RP
      gamma = HALF
      delta = 1.1E0_RP
      rho = rhobeg
      parmu = ZERO
!if (iprint >= 2) then
!print 10, RHO
!10  format(/3X, 'The initial value of RHO is', 1PE13.6, 2X, 'and PARMU is set to zero.')
!end if
      nf = 0
      sim = ZERO
      simi = ZERO
      do i = 1, n
          sim(i, i) = rho
          simi(i, i) = ONE / rho
      end do
      sim(:, n + 1) = x

      jdrop = n + 1
      ibrnch = 0

      nsav = 0
      datsav = hugenum
      almost_infinity = huge(ZERO) / TWO

40    do i = 1, n
          if (is_nan(x(i))) then
              f = x(i)
              ! Set F to NaN
              INFO = -1
              goto 600
          end if
      end do

      call calcfc(n, m, x, f, con)
      nf = nf + 1
      resmax = maxval([ZERO, -con(1:m)])
      con(m + 1) = f
      con(m + 2) = resmax

! CONSAV always containts the containt value of the current x. CON, however, will be changed during
! the calculation (see the lines above line number 220).
      consav = con

!if (nf == IPRINT - 1 .or. IPRINT == 3) then
!    print 70, nf, F, RESMAX, (X(I), I=1, IPTEM)
!70  format(/3X, 'nf =', I5, 3X, 'F =', 1PE13.6, 4X, 'MAXCV =', 1PE13.6 / 3X, 'X =', 1PE13.6, 1P4E15.6)
!    if (IPTEM < N) print 80, (X(I), I=IPTEM + 1, N)
!80  format(1PE19.6, 1P4E15.6)
!end if

! If the objective function value or the constraints contain a NaN or an infinite value, the exit.
      if (is_nan(F) .or. F > ALMOST_INFINITY) then
          info = -2
          goto 600
      end if
      if (any(is_nan(con(1:m)))) then
          resmax = sum(abs(con(1:m)))
          ! Set RESMAX to NaN
          info = -2
          goto 600
      end if
! If the objective function achieves the target value at a feasible point, then exit.
      if (f <= ftarget .and. resmax < ctol) then
          info = 1
          return
      end if

      if (nf >= maxfun .and. nf > 0) then
!    if (IPRINT >= 1) print 50
!50  format(/3X, 'Return from subroutine COBYLA because the ', 'MAXFUN limit has been reached.')
          info = 3
      end if

      if (ibrnch == 1) then
          write (10, *) '286 g440'
          goto 440
      end if

! Set the recently calculated function values in a column of DATMAT. This array has a column for
! each vertex of the current simplex, the entries of each column being the values of the constraint
! functions (if any) followed by the objective function and the greatest constraint violation at
! the vertex.
      datmat(:, jdrop) = con

      if (nf > n + 1) then
          write (10, *) '302 g130'
          goto 130
      end if

! Exchange the new vertex of the initial simplex with the optimal vertex if necessary. Then, if the
! initial simplex is not complete, pick its next vertex and calculate the function values there.
      if (jdrop <= n) then
          if (datmat(m + 1, n + 1) <= f) then
! When nf<=N, it is necessary to update X(JDROP) because next X will be calculated based on the
! current one (see the code below line number 120). The purpose of this update is to make this X
! hold the variable that has the smallest function value up to now. The next X is defined as a
! perturbation of this optimal X, which is reasonable. However, this update leads to inconsistency
! between X and [F, RESMAX]. This can cause COBYLA return inconsistent [X, F, RESMAX]. Fortunately,
! this is not a problem if NF <= N. Because, if COBYLA returns with a NF <= N, then X contains NaN
! or F = NaN or nearly Inf or the constraint contains NaN, all of which would lead to an immediate
! jump to line 600 without coming here. Therefore, if we restrict the update to only the case with
! NF <= N, ther will be no inconsistency at return.
! With the original code, if COBYLA returns with NF = N+1 (this can happen if the initial
! trust-region subproblem constantly produces too short steps so that RHO is reduced to RHOEND
! without producing any acceptable trial step; see the code below line number 380), then, when
! the code arrives at line number 600, X and [F, RESMAX] may be inconsistent. However, recall that
! the inconsistency occurs only if SIM(:, N+1) has a lower function value than X (before updated).
! Therefore, as long as the code takes SIM(:, N+1) instead of X, no harm would be done. It is the
! case likely the in the original code, because if COBYLA returns with NF=N+1, then PARMU=0, and
! hence SIM(:, N+1) will be returned.
!        x(jdrop)=sim(jdrop,n+1)
              if (nf <= n) then
                  x(jdrop) = sim(jdrop, n + 1)
              end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else
              sim(jdrop, n + 1) = x(jdrop)
              datmat(:, jdrop) = datmat(:, n + 1)
              datmat(:, n + 1) = con
              sim(jdrop, 1:jdrop) = -rho
              do k = 1, jdrop
                  simi(jdrop, k) = -sum(simi(k:jdrop, k))
              end do
          end if
      end if

! 120
      if (nf <= n) then
          jdrop = nf
          x(jdrop) = x(jdrop) + rho
          write (10, *) '367 g40'
          goto 40
      end if
130   ibrnch = 1

! Identify the optimal vertex of the current simplex.
140   jopt = n + 1
      phi = datmat(m + 1, :) + parmu * datmat(m + 2, :)
      phimin = minval(phi)
      if (phimin < phi(jopt)) then
          jopt = int(minloc(phi, dim=1), kind(jopt))
      end if
      if (parmu <= ZERO .and. minval(datmat(m + 2, :), mask=(phi <= phim&
     &in)) < datmat(m + 2, jopt)) then
! (PARMU <= ZERO) is indeed (PARMU == ZERO), and (PHI <= PHIMIN) is indeed (PHI == PHIMIN). We
! write them in this way to avoid equality comparison of real numbers.
          jopt = int(minloc(datmat(m + 2, :), mask=(phi <= phimin), dim=&
     &1), kind(jopt))
      end if

! Switch the best vertex into pole position if it is not there already, and also update SIM, SIMI
! and DATMAT.
      if (jopt <= n) then
          datmat(:, [jopt, n + 1]) = datmat(:, [n + 1, jopt])
          ! Exchange DATMAT(:, JOPT) AND DATMAT(:, N+1)
          sim(:, n + 1) = sim(:, n + 1) + sim(:, jopt)
          simjopt = sim(:, jopt)
          sim(:, jopt) = ZERO
          sim(:, 1:n) = sim(:, 1:n) - spread(simjopt, dim=2, ncopies=n)
! The above update is equivalent to multiply SIM(:, 1:N) from the right side by a matrix whose
! JOPT-th row is [-1, -1, ..., -1], while all the other rows are the same as those of the
! identity matrix. It is easy to check that the inverse of this matrix is itself. Therefore,
! SIMI should be updated by a multiplication with this matrix (i.e., its inverse) from the left
! side, as is done in the following line. The JOPT-th row of the updated SIMI is minus the sum
! of all rows of the original SIMI, whereas all the other rows remain unchanged.
          simi(jopt, :) = -sum(simi, dim=1)
      end if

      if (info == 3) then
          write (10, *) '420 g600'
          goto 600
      end if

! Make an error return if SIMI is a poor approximation to the inverse of SIM(:, 1:N).
      error = matprod(simi, sim(:, 1:n))
      do i = 1, n
          error(i, i) = -ONE + error(i, i)
      end do
      if (any(is_nan(error)) .or. maxval(abs(error)) > TENTH) then
!    if (IPRINT >= 1) print 210
!210 format(/3X, 'Return from subroutine COBYLA because ', 'rounding errors are becoming damaging.')
          info = 7
          write (10, *) '457 g600'
          goto 600
      end if

! Calculate the coefficients of the linear approximations to the objective and constraint functions,
! placing minus the objective function gradient after the constraint gradients in the array A.
! When __USE_INTRINSIC_ALGEBRA__ = 1, the following code may not produce the same result as
! Powell's, because the intrinsic MATMUL behaves differently from a naive triple loop in
! finite-precision arithmetic.
! 220
      con = -datmat(:, n + 1)
      ! Why put a negative sign?????
! Is it more reasonable to save A transpose instead of A? Better name for A?
      A = transpose(matprod(datmat(1:m + 1, 1:n) - spread(datmat(1:m + 1&
     &, n + 1), dim=2, ncopies=n), simi))
      A(:, m + 1) = -A(:, m + 1)

! Calculate the values of sigma and eta, and set IFLAG=0 if the current simplex is not acceptable.
      parsig = alpha * rho
      pareta = beta * rho
! For J = 1, 2, ..., n, VSIG(J) is The Euclidean distance from vertex J to the opposite face of
! the current simplex. But what about vertex N+1?
      vsig = ONE / sqrt(sum(simi**2, dim=2))
      veta = sqrt(sum(sim(:, 1:n)**2, dim=1))
!---------------------------------------------------------------------------------------!
!improve_geo = any(vsig < parsig) .or. any(veta > pareta) .or. any(is_nan([vsig, veta]))
!---------------------------------------------------------------------------------------!
      improve_geo = any(vsig < parsig) .or. any(veta > pareta)

      if (IBRNCH == 1 .or. .not. improve_geo) then
          write (10, *) '515 g370'
          goto 370
      end if

! If a new vertex is needed to improve acceptability, then decide which vertex to drop from the simplex.
      if (maxval(veta) > pareta) then
          jdrop = maxloc(veta, dim=1)
      elseif (minval(vsig) < pareta) then
          jdrop = minloc(vsig, dim=1)
      else
          jdrop = 0_IK
      end if

! If VETA or VSIG become NaN due to rounding errors, JDROP will end up being 0. If we continue, then
! a Segmentation Fault will happen because we will access SIM(:, JDROP) and VSIG(JDROP).
      if (jdrop == 0) then
!    if (IPRINT >= 1) print 286
!286 format(/3X, 'Return from subroutine COBYLA because ', 'rounding errors are becoming damaging.')
          INFO = 7
          goto 600
      end if

! Save the information of the JOPT-th vertex in XSAV and DATSAV.
      call savex(sim(:, n + 1) + sim(:, jdrop), datmat(:, jdrop), xsav, &
     &datsav, n, m, nsav, nsmax, ctol)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the step to the new vertex and its sign.
!
      TEMP = GAMMA * RHO * VSIG(jdrop)
      do I = 1, N
          DX(I) = TEMP * SIMI(jdrop, I)
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      CVMAXP=0.0
!      CVMAXM=0.0
      CVMAXP = ZERO
      CVMAXM = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do K = 1, M + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SUMM=0.0
          SUMM = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              SUMM = SUMM + A(I, K) * DX(I)
          end do
          if (K < M + 1) then
              TEMP = DATMAT(K, N + 1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          CVMAXP=AMAX1(CVMAXP,-SUMM-TEMP)
!          CVMAXM=AMAX1(CVMAXM,SUMM-TEMP)
              CVMAXP = DMAX1(CVMAXP, -SUMM - TEMP)
              CVMAXM = DMAX1(CVMAXM, SUMM - TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          end if
      end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      DXSIGN=1.0
!      IF (PARMU*(CVMAXP-CVMAXM) .GT. SUMM+SUMM) DXSIGN=-1.0
      DXSIGN = ONE
      if (PARMU * (CVMAXP - CVMAXM) > SUMM + SUMM) DXSIGN = -ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Update the elements of SIM and SIMI, and set the next X.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
      TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          DX(I) = DXSIGN * DX(I)
          SIM(I, jdrop) = DX(I)
          TEMP = TEMP + SIMI(jdrop, I) * DX(I)
      end do
      do I = 1, N
          SIMI(jdrop, I) = SIMI(jdrop, I) / TEMP
      end do
      do J = 1, N
          if (J /= jdrop) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=0.0
              TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do I = 1, N
                  TEMP = TEMP + SIMI(J, I) * DX(I)
              end do
              do I = 1, N
                  SIMI(J, I) = SIMI(J, I) - TEMP * SIMI(jdrop, I)
              end do
          end if
          X(J) = SIM(J, N + 1) + DX(J)
      end do
      write (10, *) '624 g40'
      goto 40
!
!     Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO.
!
370   IZ = 1
      IZDOTA = IZ + N * N
      IVMC = IZDOTA + N
      ISDIRN = IVMC + M + 1
      IDXNEW = ISDIRN + N
      IVMD = IDXNEW + N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2019-08-29: For ill-conditioned problems, NaN may occur in the
! models. In such a case, we terminate the code. Otherwise, the behavior
! of TRSTLP is not predictable, and Segmentation Fault or infinite
! cycling may happen. This is because any equality/inequality comparison
! involving NaN returns FALSE, which can lead to unintended behavior of
! the code, including uninitialized indices.
      do J = 1, N
          do I = 1, N
              if (is_nan(SIMI(I, J))) then
                  if (IPRINT >= 1) print 376
376   format(/3X, 'Return from subroutine COBYLA because ', 'rounding er&
     &rors are becoming damaging.')
                  INFO = 7
                  write (10, *) '648 g600'
                  goto 600
              end if
          end do
      end do
      do J = 1, M + 1
          do I = 1, N
              if (is_nan(A(I, J))) then
                  INFO = -3
                  write (10, *) '657 g600'
                  goto 600
              end if
          end do
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write (10, *) nf, 'bt'
!call TRSTLP(N, M, A, CON, RHO, DX, IFULL, IACT)!, W(IZ), W(IZDOTA), W(IVMC), W(ISDIRN), W(IDXNEW))!, W(IVMD:))
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!call TRSTLP(N, M, A, CON, RHO, DX, IFULL, IACT, W(IVMD:))
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


      if (nf == 8) then
          write (12, *), A(1:n, 1:m + 1), con(1:m + 2)
          close (12)
      end if
      call TRSTLP(N, M, A, CON, RHO, DX, IFULL, IACT)
      write (10, *) nf, 'at'
      if (IFULL == 0) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=0.0
          TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              TEMP = TEMP + DX(I)**2
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .LT. 0.25*RHO*RHO) THEN
          if (TEMP < 0.25D0 * RHO * RHO) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              IBRNCH = 1
              write (10, *) '657 g550'
              goto 550
          end if
      end if
!
!     Predict the change to F and the new maximum constraint violation if the
!     variables are altered from x(0) to x(0)+DX.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RESNEW=0.0
!      CON(MP)=0.0
      RESNEW = ZERO
      CON(M + 1) = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do K = 1, M + 1
          SUMM = CON(K)
          do I = 1, N
              SUMM = SUMM - A(I, K) * DX(I)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (K .LT. MP) RESNEW=AMAX1(RESNEW,SUMM)
          if (K < M + 1) RESNEW = DMAX1(RESNEW, SUMM)
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
      BARMU = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PREREC = DATMAT(m + 2, N + 1) - RESNEW
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (PREREC .GT. 0.0) BARMU=SUMM/PREREC
!      IF (PARMU .LT. 1.5*BARMU) THEN
!          PARMU=2.0*BARMU
      if (PREREC > ZERO) BARMU = SUMM / PREREC
      if (PARMU < 1.5D0 * BARMU) then
          PARMU = TWO * BARMU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (IPRINT >= 2) print 410, PARMU
410   format(/3X, 'Increase in PARMU to', 1PE13.6)
          PHITMP = DATMAT(M + 1, N + 1) + PARMU * DATMAT(m + 2, N + 1)
          do J = 1, N
              TEMP = DATMAT(M + 1, J) + PARMU * DATMAT(m + 2, J)
              if (TEMP < PHITMP) then
                  write (10, *) '706 g140'
                  goto 140
              end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .EQ. PHITMP .AND. PARMU .EQ. 0.0) THEN
              if (abs(TEMP - PHITMP) <= ZERO .and. PARMU <= ZERO) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  if (DATMAT(m + 2, J) < DATMAT(m + 2, N + 1)) then
                      write (10, *) '714 g140'
                      goto 140
                  end if
              end if
          end do
      end if
      PREREM = PARMU * PREREC - SUMM
!
!     Calculate the constraint and objective functions at x(*). Then find the
!     actual reduction in the merit function.
!
      do I = 1, N
          X(I) = SIM(I, N + 1) + DX(I)
      end do
      IBRNCH = 1
      write (10, *) '729 g40'
      goto 40
440   VMOLD = DATMAT(M + 1, N + 1) + PARMU * DATMAT(m + 2, N + 1)
      VMNEW = F + PARMU * RESMAX
      TRURED = VMOLD - VMNEW
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (PARMU .EQ. 0.0 .AND. F .EQ. DATMAT(MP,N+1)) THEN
      if (PARMU <= ZERO .and. abs(F - DATMAT(M + 1, N + 1)) <= ZERO) the&
     &n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          PREREM = PREREC
          TRURED = DATMAT(m + 2, N + 1) - RESMAX
      end if
!
!     Begin the operations that decide whether x(*) should replace one of the
!     vertices of the current simplex, the change being mandatory if TRURED is
!     positive. Firstly, jdrop is set to the index of the vertex that is to be
!     replaced.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RATIO=0.0
!      IF (TRURED .LE. 0.0) RATIO=1.0
      RATIO = ZERO
      if (TRURED <= ZERO) RATIO = ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      jdrop = 0
      do J = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
          TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do I = 1, N
              TEMP = TEMP + SIMI(J, I) * DX(I)
          end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=ABS(TEMP)
          TEMP = DABS(TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (TEMP > RATIO) then
              jdrop = J
              RATIO = TEMP
          end if
          SIGBAR(J) = TEMP * VSIG(J)
      end do
!
!     Calculate the value of ell.
!
      EDGMAX = DELTA * RHO
      L = 0
      do J = 1, N
          if (SIGBAR(J) >= PARSIG .or. SIGBAR(J) >= VSIG(J)) then
              TEMP = VETA(J)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TRURED .GT. 0.0) THEN
!              TEMP=0.0
              if (TRURED > ZERO) then
                  TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  do I = 1, N
                      TEMP = TEMP + (DX(I) - SIM(I, J))**2
                  end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              TEMP=SQRT(TEMP)
                  TEMP = DSQRT(TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              end if
              if (TEMP > EDGMAX) then
                  L = J
                  EDGMAX = TEMP
              end if
          end if
      end do
      if (L > 0) jdrop = L
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820:
! When jdrop=0, the algorithm decides not to include the trust-region
! trial point X into the simplex, because X is not good enough according
! to the merit function PHI = F + PARMU*RESMAX. In this case, X will
! simply be discarded in the original code. However, this decision
! depends on the value of PARMU. When PARMU is updated later, the
! discarded X might turn out better, sometimes even better than
! SIM(:, N+1), which is supposed to be the best point in the simplex.
! For this reason, we save the to-be-discarded X in XSAV and compare
! them with SIM(:, N+1) right before exiting. If a vector in XSAV turns
! out better than SIM(:, N+1), we replace SIM(:, N+1) by this vector
!
! When jdrop > 0, SIM(:, jdrop) will be removed from the simplex
! according to PHI with the current PARMU. Similar to X, SIM(:, jdrop)
! may turn out better when PARMU is updated. Therefore, XSAV also takes
! SIM(:, jdrop) into account.
!
! We save at most NSMAX to-be-discarded X.
!
      if (jdrop == 0) then
          call savex(x, consav, xsav, datsav, n, m, nsav, nsmax, ctol)
      else
          call savex(sim(:, n + 1) + sim(:, jdrop), datmat(:, jdrop), xs&
     &av, datsav, n, m, nsav, nsmax, ctol)
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (jdrop == 0) then
          write (10, *) '841 g550'
          goto 550
      end if
!
!     Revise the simplex by updating the elements of SIM, SIMI and DATMAT.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
      TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do I = 1, N
          SIM(I, jdrop) = DX(I)
          TEMP = TEMP + SIMI(jdrop, I) * DX(I)
      end do
      do I = 1, N
          SIMI(jdrop, I) = SIMI(jdrop, I) / TEMP
      end do
      do J = 1, N
          if (J /= jdrop) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=0.0
              TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do I = 1, N
                  TEMP = TEMP + SIMI(J, I) * DX(I)
              end do
              do I = 1, N
                  SIMI(J, I) = SIMI(J, I) - TEMP * SIMI(jdrop, I)
              end do
          end if
      end do
      do K = 1, m + 2
          DATMAT(K, jdrop) = CON(K)
      end do
!
!     Branch back for further iterations with the current RHO.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (TRURED .GT. 0.0 .AND. TRURED .GE. 0.1*PREREM) GOTO 140
      if (TRURED > ZERO .and. TRURED >= TENTH * PREREM) then
          write (10, *) '881 g140'
          goto 140
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
550   if (improve_geo) then
          IBRNCH = 0
          write (10, *), '886 g140'
          goto 140
      end if
!
!     Otherwise reduce RHO if it is not at its least value and reset PARMU.
!
      if (RHO > RHOEND) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          RHO=0.5*RHO
!          IF (RHO .LE. 1.5*RHOEND) RHO=RHOEND
!          IF (PARMU .GT. 0.0) THEN
!              DENOM=0.0
          RHO = 0.5D0 * RHO
          if (RHO <= 1.5D0 * RHOEND) RHO = RHOEND
          if (PARMU > ZERO) then
              DENOM = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do K = 1, M + 1
                  CMIN = DATMAT(K, N + 1)
                  CMAX = CMIN
                  do I = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              CMIN=AMIN1(CMIN,DATMAT(K,I))
!  560         CMAX=AMAX1(CMAX,DATMAT(K,I))
!              IF (K .LE. M .AND. CMIN .LT. 0.5*CMAX) THEN
!                  TEMP=AMAX1(CMAX,0.0)-CMIN
!                  IF (DENOM .LE. 0.0) THEN
                      CMIN = DMIN1(CMIN, DATMAT(K, I))
                      CMAX = DMAX1(CMAX, DATMAT(K, I))
                  end do
                  if (K <= M .and. CMIN < 0.5D0 * CMAX) then
                      TEMP = DMAX1(CMAX, ZERO) - CMIN
                      if (DENOM <= ZERO) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          DENOM = TEMP
                      else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                      DENOM=AMIN1(DENOM,TEMP)
                          DENOM = DMIN1(DENOM, TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      end if
                  end if
              end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              IF (DENOM .EQ. 0.0) THEN
!                  PARMU=0.0
              if (abs(DENOM) <= ZERO) then
              ! DENOM <= ZERO???  Is it nonnegative?
                  PARMU = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              else if (CMAX - CMIN < PARMU * DENOM) then
                  PARMU = (CMAX - CMIN) / DENOM
              end if
          end if
          if (IPRINT >= 2) print 580, RHO, PARMU
580   format(/3X, 'Reduction in RHO to', 1PE13.6, ' and PARMU =', 1PE13.&
     &6)
!    if (IPRINT == 2) then
!        print 70, nf, DATMAT(M + 1, N + 1), DATMAT(m + 2, N + 1), (SIM(I, N + 1), I=1, IPTEM)
!        if (IPTEM < N) print 80, (X(I), I=IPTEM + 1, N)
!    end if
          write (10, *), '946 g140'
          goto 140
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      else
          INFO = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if
!
!     Return the best calculated values of the variables.
!
      if (IPRINT >= 1) print 590
590   format(/3X, 'Normal return from subroutine COBYLA')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (IFULL .EQ. 1) GOTO 620
!  600 DO 610 I=1,N
!  610 X(I)=SIM(I,N+1)
!      F=DATMAT(MP,N+1)
!      RESMAX=DATMAT(m+2,N+1)
!
!      Zaikun 01-06-2019:
!      Why go to 620 directly without setting X and F? This seems
!      INCORRECT, because it may lead to a return with X and F that
!      are not the best available.
!      The following code defines X as an "optimal" one among the
!      vectors:
!      DATSAV(:, 1:NSAV), (when NSAV >= 1),
!      SIM(:, N+1), and SIM(:, 1:min(N, nf-2)) (when nf>=2).
!      Here, X being "optimal" means
!      1. the constraint violation of X is at most RESREF
!      2. no other vector is better than X according to ISBETTER with
!      the current PARMU.
!
!      Note:
!      0. The last evaluated X and its function/constraint information
!      are saved in [X, CONSAV, F, RESMAX].
!      1. When nf=1, SIM and DATMAT have not been initialized yet.
!      2. When 2<=nf<=N+1, the first evaluated X are saved in
!      SIM(:, N+1), its function/constraint in DATMAT(:, N+1), while the
!      other X are saved in SIM(:, nf-1), its function/constraint
!      in DATMAT(:, nf-1). However, when the code arrives at line 600,
!      [X, CON, F, RESMAX] may have not been saved into SIM(:, nf-1)
!      and DATMAT(:, nf-1) yet. That is why we check SIM up to
!      nf-2 instead of nf-1.
600   con = consav
      parmu = max(parmu, 1.0E2_RP)
      if (nf >= 2) then
      ! See the comments above for why NF > 2.
          if (isbetter(f, resmax, datmat(m + 1, n + 1), datmat(m + 2, n &
     &+ 1), parmu, ctol)) then
              x = sim(:, n + 1)
              f = datmat(m + 1, n + 1)
              resmax = datmat(m + 2, n + 1)
              con = datmat(:, n + 1)
          end if
          resref = resmax
          ! We only consider points with constraint violation no more than RESREF.
          if (is_nan(resref)) then
              resref = HUGENUM
          end if
          do j = 1, min(n, nf - 2)
          ! See the comments above for why to check these J.
              if (datmat(m + 2, j) <= resref .and. isbetter(f, resmax, d&
     &atmat(m + 1, j), datmat(m + 2, j), parmu, ctol)) then
                  x = sim(:, j) + sim(:, n + 1)
                  f = datmat(m + 1, j)
                  resmax = datmat(m + 2, j)
                  con = datmat(:, j)
              end if
          end do
      end if

!jmax = min(n, nf - 2)
!if (minval(datmat(m + 1, 1:jmax) + parmu * datmat(m + 2, 1:jmax), mask=(datmat(m + 2, 1:jmax) <= resref)) < f + parmu * resmax) then
!    j = int(minloc(datmat(m + 1, 1:jmax) + parmu * datmat(m + 2, 1:jmax), mask=(datmat(m + 2, 1:jmax) <= resref), dim=1), kind(j))
!    x = sim(:, j) + sim(:, n + 1)
!    f = datmat(m + 1, j)
!    resmax = datmat(m + 2, j)
!    con = datmat(:, j)
!end if

!if (minval(datsav(m + 1, 1:nsav) + parmu * datsav(m + 2, 1:nsav), mask=(datsav(m + 2, 1:nsav) <= resref)) < f + parmu * resmax) then
!    j = int(minloc(datsav(m + 1, 1:nsav) + parmu * datsav(m + 2, 1:nsav), mask=(datsav(m + 2, 1:nsav) <= resref), dim=1), kind(j))
!    x = xsav(:, j)
!    f = datsav(m + 1, j)
!    resmax = datsav(m + 2, j)
!    con = datsav(:, j)
!end if
      do j = 1, nsav
          if (datsav(m + 2, j) <= resref .and. isbetter(f, resmax, datsa&
     &v(m + 1, j), datsav(m + 2, j), parmu, ctol)) then
              x = xsav(:, j)
              f = datsav(m + 1, j)
              resmax = datsav(m + 2, j)
              con = datsav(:, j)
          end if
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!620 if (IPRINT >= 1) then
!    print 70, nf, F, RESMAX, (X(I), I=1, IPTEM)
!    if (IPTEM < N) print 80, (X(I), I=IPTEM + 1, N)
!end if

      return
      end subroutine cobylb

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820: See the comments below line number 480
      subroutine SAVEX(XDROP, DATDROP, XSAV, DATSAV, N, M, NSAV, NSMAX, &
     &CTOL)
! This subroutine saves XDROP in XSAV and DATDROP in DATSAV, unless
! XDROP is dominated by a vector in XSAV(:, 1:NSAV). If XDROP dominates
! some vectors in XSAV(:, 1:NSAV), then these vectors will be removed.
! If XDROP does not dominate any of XSAV(:, 1:NSAV) but NSAV=NSMAX,
! then we remove XSAV(:,1), which is the oldest vector in XSAV(:, 1:NSAV).
!
! When COBYLA calls this subroutine, XDROP is a vector to be "dropped",
! and  DATDROP contains its function/constraint inforation (in particular,
! DATDROP(MP) = F(XDROP), and DATDROP(m+2) = RESMAX(X)). XSAV and DATSAV
! save at most NSMAX vectors "dropped" by COBYLB and their function/constraint
! information. Only XSAV(:, 1:NSAV) and DATSAV(:, 1:NSAV) contains such
! vectors, while XSAV(:, NSAV+1:NSMAX) and DATSAV(:, NSAV+1:NSMAX) are
! not initialized yet.
!
! Note: X dominates Y if and only if the function/constraint of X is
! better than the function/constraint of Y accoring to the ISBETTER
! subroutine with PARMU = -ONE. Indeed, PARMU can be any negative
! number. This is because, due to the implementation of ISBETTER,
! X dominates Y (i.e., X is better than Y with PARMU < 0)
! ==> X is better than Y with any PARMU >= 0,
! ==> X is better than Y regardless of PARMU.
! For this reason, it is sufficient to save all the "dropped" X that are not
! dominated by any vectors in XSAV (as we do in this subroutine),
! because the other X are always worse than some vector in XSAV.
!

! Generic modules
      use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, HUGENU&
     &M, DEBUGGING, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : calquad, inprod

! Solver-specific modules
!use isbetter_mod, only : isbetter

      implicit none

! Inputs
      integer(IK), intent(IN) :: m
      integer(IK), intent(IN) :: n
      integer(IK), intent(IN) :: nsmax
      real(RP), intent(IN) :: ctol
      real(RP), intent(IN) :: datdrop(:)
      ! m+2
      real(RP), intent(IN) :: xdrop(:)
      ! n

! In-outputs
      integer(IK), intent(INOUT) :: Nsav
      real(RP), intent(INOUT) :: Datsav(:, :)
      ! (M+2, NSMAX)
      real(RP), intent(INOUT) :: Xsav(:, :)
      ! (N, NSMAX)

! Local variables
      integer(IK) :: i
      integer(IK) :: j
      integer(IK) :: k
      integer(IK) :: l
      integer(IK) :: mp
      integer(IK) :: mpp
      integer(IK) :: nremove
      integer(IK) :: iremove(nsmax)
      real(RP) :: parmu
      logical :: better

      if (NSMAX <= 0) return
      ! Do nothing if NSMAX=0

      PARMU = -ONE
      ! See the comments above for why PARMU = -1

! IREMOVE: indices of vectors to remove from XSAV
! NREMOVE: number of vectors to remove from XSAV
      do J = 1, NSMAX
! It is not enough to initialize IREMOVE(1:NSAV), because NSAV may be
! incremented by 1 latter, and then IREMOVE(NSAV+1) will be accessed.
          IREMOVE(J) = -1
      end do
      NREMOVE = 0
      do I = 1, NSAV
! If XDROP is dominated by XSAV(:, I), then return immediately,
! because XDROP should not be inluded into XSAV.
          if (isbetter(datdrop(m + 1), datdrop(m + 2), datsav(m + 1, i),&
     & datsav(m + 2, i), parmu, ctol)) return
! If XDROP dominates XSAV(:, I), then increment NREMOVE by 1 and save
! I as IREMOVE(NREMOVE).
          if (isbetter(datsav(m + 1, i), datsav(m + 2, i), datdrop(m + 1&
     &), datdrop(m + 2), parmu, ctol)) then
              NREMOVE = NREMOVE + 1
              IREMOVE(NREMOVE) = I
          end if
      end do

! The code did not return and NREMOVE=0 (no vector to remove from XSAV).
! If NSAV=NSMAX, then remove XSAV(:, 1); otherwise, increment NSAV by
! 1 and then "remove" XSAV(:, NSAV) (though there is no vector saved there)
      if (NREMOVE == 0) then
          if (NSAV == NSMAX) then
              IREMOVE(1) = 1
          else
              NSAV = NSAV + 1
              IREMOVE(1) = NSAV
          end if
          NREMOVE = 1
      end if

! Remove from XSAV the vectors indexed by IREMOVE
      J = 1
      ! Index of IREMOVE
      K = 1
      ! Index of the new XSAV
      do I = 1, NSAV
      ! Index of the old XSAV
          if (I == IREMOVE(J)) then
              J = J + 1
              ! Do nothing but incrementing J by 1
          else
          ! Set XSAV(:, K) = XSAV(:, I)
              do L = 1, N
                  XSAV(L, K) = XSAV(L, I)
              end do
              do L = 1, m + 2
                  DATSAV(L, K) = DATSAV(L, I)
              end do
              K = K + 1
              ! Increment K by 1
          end if
      end do

! Set the number of vectors in the new XSAV
      NSAV = NSAV - NREMOVE + 1

! Save XDROP in XSAV(:, NSAV) (with NSAV updated as above)
      if (NSAV >= 1 .and. NSAV <= NSMAX) then
! This inequlity is not guaranteed if NSMAX=0, where NSAV will
! be 0 and hence a Segmentation Fault when accessing
! XSAV(:, NSAV). Although we return immediately if NSMAX=0,
! we still check this inequlity to be safe.
          do L = 1, N
              XSAV(L, NSAV) = XDROP(L)
          end do
          do L = 1, m + 2
              DATSAV(L, NSAV) = DATDROP(L)
          end do
      end if

      end subroutine SAVEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820:
      function isbetter(f0, conv0, f, conv, parmu, ctol) result(better)
! This subroutine compares whether (F, CONV) is (strictly) better than (F0, CONV0) in the sense of
! decreasing the merit function PHI = F + PARMU*CONV. It takes care of the cases where some of these
! values are NaN or Inf. At return, BETTER=true iff (F, CONV) is better than (F0, CONV0).
! N.B.:
! 1. A = Inf if and only if A > HUGENUM (defined below);
! 2. A = NaN if and only if A /= A;
! 3. If A = NaN, then any comparison (except /=) with another number B (can be Inf or NaN) returns false.

! Generic modules
      use consts_mod, only : RP, IK, ZERO, HALF, TENTH, HUGENUM, DEBUGGI&
     &NG, SRNLEN
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf
      use debug_mod, only : errstop
      use output_mod, only : retmssg, rhomssg, fmssg
      use lina_mod, only : calquad, inprod
      implicit none

      real(RP), intent(IN) :: f0
      real(RP), intent(IN) :: conv0
      real(RP), intent(IN) :: f
      real(RP), intent(IN) :: conv
      real(RP), intent(IN) :: parmu
      real(RP), intent(IN) :: ctol
      logical :: better

! Local variables
      logical :: f0infnan
      logical :: finfnan
      logical :: fle
      logical :: flt
      logical :: c0infnan
      logical :: cinfnan
      logical :: cle
      logical :: clt

      better = .false.

! As values of F0, CONV0, F, and CONV, we regard Inf and NaN being equivalent values (they are equally bad).
      f0infnan = is_nan(f0) .or. (f0 > HUGENUM)
      ! F0 = Inf or NaN?
      c0infnan = is_nan(conv0) .or. (conv0 > HUGENUM)
      ! CONV0 = Inf or NaN?
      finfnan = is_nan(f) .or. (f > HUGENUM)
      ! F = Inf or NaN?
      cinfnan = is_nan(conv) .or. (conv > HUGENUM)
      ! CONV = Inf or NaN?

! If PARMU >= 0 and F + PARMU*CONV < F0 + PARMU*CONV0 and CONV < CTOL (feasible), then (F, CONV) is
! better than (F0, CONV0).
! Note that we should not set BETTER=FALSE even if this inequality does not hold, as either side may be NaN.
      if (parmu >= zero .and. f + parmu * conv < f0 + parmu * conv0 .and&
     &. conv < ctol) then
          better = .true.
      end if

! If CONV < CTOL and F is not Inf or NaN while (CONV0 < CTOL) is false (may be because CONV0 is NaN),
! then (F, CONV) is better than (F0, CONV0). We prefer feasible points (i.e., constraint violation
! is less than CTOL) to infeasible ones.
      if (conv < ctol .and. .not. (conv0 < ctol) .and. .not. finfnan) th&
     &en
          better = .true.
      end if

! If F0 or CONV0 is Inf/NaN while neither F nor CONV is Inf/NaN, then (F, CONV) is better than (F0, CONV0).
      if ((f0infnan .or. c0infnan) .and. .not. (finfnan .or. cinfnan)) t&
     &hen
          better = .true.
      end if

      flt = (f0infnan .and. (.not. finfnan)) .or. (f < f0)
      ! F < F0?
      fle = (f0infnan .and. finfnan) .or. (f <= f0) .or. flt
      ! F <= F0?
      clt = (c0infnan .and. (.not. cinfnan)) .or. (conv < conv0)
      ! CONV < CONV0?
      cle = (c0infnan .and. cinfnan) .or. (conv <= conv0) .or. clt
      ! CONV <= CONV0?

! If (F < F0 and CONV <= CONV0) or (F <= F0 and CONV < CONV0) in the sense defined above, the
! (F, CONV) is better than (F0, CONV0).
      if ((flt .and. cle) .or. (fle .and. clt)) then
          better = .true.
      end if

! If one of F and CONV is -Inf and the other is not Inf/Nan, while neither F0 nor CONV0 is -Inf,
! then the (F, CONV) is better than (F0, CONV0).
      if ((f < -HUGENUM) .and. (.not. cinfnan) .and. (f0 >= -HUGENUM) .a&
     &nd. (conv0 >= -HUGENUM)) then
          better = .true.
      end if
      if ((conv < -HUGENUM) .and. (.not. finfnan) .and. (f0 >= -HUGENUM)&
     & .and. (conv0 >= -HUGENUM)) then
          better = .true.
      end if
      end function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module cobylb_mod