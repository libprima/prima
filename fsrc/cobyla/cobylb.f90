module cobylb_mod

contains

subroutine cobylb(n, m, x, rhobeg, rhoend, iprint, maxfun, con, f, info, ftarget, resmax)

! Generic modules
use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, HUGENUM, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : calquad, inprod

! Solver-specific modules
!use savex_mod, only : savex
!use isbetter_mod, only : isbetter
use trustregion_mod, only : trstlp

implicit none

! Inputs
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: m
integer(IK), intent(in) :: maxfun
integer(IK), intent(in) :: n
real(RP), intent(in) :: ftarget
real(RP), intent(in) :: rhobeg
real(RP), intent(in) :: rhoend

! In-outputs
real(RP), intent(inout) :: con(:) ! m+2. Bad name; should be confr
real(RP), intent(inout) :: x(:)  ! n

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
integer(IK) :: iflag
integer(IK) :: ifull
integer(IK) :: iptem
integer(IK) :: iptemp
integer(IK) :: isdirn
integer(IK) :: ivmc
integer(IK) :: ivmd
integer(IK) :: iz
integer(IK) :: izdota
integer(IK) :: j
integer(IK) :: jdrop
integer(IK) :: k
integer(IK) :: l
integer(IK) :: mp
integer(IK) :: mpp
integer(IK) :: nbest
integer(IK) :: nfvals
integer(IK) :: np
integer(IK) :: nsav
logical :: better
real(RP) :: A(size(x), m + 1)  ! (n, )
real(RP) :: almost_infinity
real(RP) :: alpha
real(RP) :: barmu
real(RP) :: beta
real(RP) :: cmax
real(RP) :: cmin
real(RP) :: consav(size(con))
real(RP) :: csum
real(RP) :: cvmaxm
real(RP) :: cvmaxp
real(RP) :: datdrop(size(con))
real(RP) :: datmat(m + 2, n + 1)  !(mpp, )
real(RP) :: datsav(size(con), nsmax)
real(RP) :: delta
real(RP) :: denom
real(RP) :: dx(size(x))
real(RP) :: dxsign
real(RP) :: edgmax
real(RP) :: error
real(RP) :: gamma
real(RP) :: pareta
real(RP) :: parmu
real(RP) :: parsig
real(RP) :: phi
real(RP) :: phimin
real(RP) :: prerec
real(RP) :: prerem
real(RP) :: ratio
real(RP) :: resmax
real(RP) :: resnew
real(RP) :: resref
real(RP) :: rho
real(RP) :: sigbar(size(x))
real(RP) :: sim(size(x), size(x) + 1)  ! (n, )
real(RP) :: simi(size(x), size(x) + 1)  ! (n, )
real(RP) :: sum
real(RP) :: temp
real(RP) :: tempa
real(RP) :: trured
real(RP) :: veta(size(x))
real(RP) :: vmnew
real(RP) :: vmold
real(RP) :: vsig(size(x))
real(RP) :: w(N * (3 * N + 2 * M + 11) + 4 * M + 6)
real(RP) :: weta
real(RP) :: wsig
real(RP) :: xdrop(size(x))
real(RP) :: xsav(size(x), NSMAX)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set the initial values of some parameters. The last column of SIM holds
!     the optimal vertex of the current simplex, and the preceding N columns
!     hold the displacements from the optimal vertex to the other vertices.
!     Further, SIMI holds the inverse of the matrix that is contained in the
!     first N columns of SIM.
!
INFO = 2147483647
IPTEM = MIN0(N, 5)
IPTEMP = IPTEM + 1
NP = N + 1
MP = M + 1
MPP = M + 2
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ALPHA=0.25
!      BETA=2.1
!      GAMMA=0.5
!      DELTA=1.1
ALPHA = 0.25D0
BETA = 2.1D0
GAMMA = 0.5D0
DELTA = 1.1D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RHO = RHOBEG
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      PARMU=0.0
PARMU = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (IPRINT >= 2) print 10, RHO
10 format(/3X, 'The initial value of RHO is', 1PE13.6, 2X, 'and PARMU is set to zero.')
NFVALS = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=1.0/RHO
TEMP = ONE / RHO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do I = 1, N
    SIM(I, NP) = X(I)
    do J = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SIM(I,J)=0.0
!   20 SIMI(I,J)=0.0
        SIM(I, J) = ZERO
        SIMI(I, J) = ZERO
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SIM(I, I) = RHO
    SIMI(I, I) = TEMP
end do
JDROP = NP
IBRNCH = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
NSAV = 0
DATSAV = HUGENUM
ALMOST_INFINITY = huge(ZERO) / TWO
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
40 do I = 1, N
    if (is_nan(X(I))) then
        F = X(I) ! Set F to NaN
        INFO = -1
        write (10, *), '200 g600'
        goto 600
    end if
end do

call CALCFC(N, M, X, F, CON)
NFVALS = NFVALS + 1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      RESMAX=0.0
RESMAX = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (M > 0) then
    do K = 1, M
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   60     RESMAX=AMAX1(RESMAX,-CON(K))
        RESMAX = DMAX1(RESMAX, -CON(K))
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if


if (NFVALS == IPRINT - 1 .or. IPRINT == 3) then
    print 70, NFVALS, F, RESMAX, (X(I), I=1, IPTEM)
70  format(/3X, 'NFVALS =', I5, 3X, 'F =', 1PE13.6, 4X, 'MAXCV =', 1PE13.6 / 3X, 'X =', 1PE13.6, 1P4E15.6)
    if (IPTEM < N) print 80, (X(I), I=IPTEMP, N)
80  format(1PE19.6, 1P4E15.6)
end if
CON(MP) = F
CON(MPP) = RESMAX

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! By Zaikun 20190819:
! CONSAV always containts the containt value of the current x.
! CON, however, will be changed during the calculation (see the lines
! above line number 220).
do K = 1, MPP
    CONSAV(K) = CON(K)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom/Zaikun (on 04-06-2019/10-06-2019):
!     CSUM containts the sum of the absolute value of the constraints to
!     check whether it contains a NaN value.
CSUM = ZERO
do K = 1, M
    CSUM = CSUM + DABS(CON(K))
end do
if (is_nan(CSUM)) then
    RESMAX = CSUM ! Set RESMAX to NaN
    INFO = -2
    write (10, *) '252 g600'
    goto 600
end if
!     If the objective function value or the constraints contain a NaN or an
!     infinite value, the algorithm stops.
if (is_nan(F) .or. F > ALMOST_INFINITY) then
    INFO = -2
    write (10, *) '258 g600'
    goto 600
end if
!     If the objective function achieves the target value at a feasible
!     point, then exit.
!      IF (F .LE. FTARGET .AND. RESMAX .LE. ZERO) THEN
if (F <= FTARGET .and. RESMAX < CTOL) then
!         The feasibility is guarantee because RESMAX .LE. CTOL
    INFO = 1
    write (10, *) '268 g620'
    goto 620
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Zaikun (on 06-06-2019)
!     The following code was placed before "CALL CALCFC (N,M,X,F,CON)".
!     This led to a bug, because F may not equal F(X) if the subroutine
!     exits due to NFVALS .GE. MAXFUN (X is updated but F is not evaluated
!     at X). Similar thing can be said about RESMAX.
if (NFVALS >= MAXFUN .and. NFVALS > 0) then
    if (IPRINT >= 1) print 50
50  format(/3X, 'Return from subroutine COBYLA because the ', 'MAXFUN limit has been reached.')
    INFO = 3
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (IBRNCH == 1) then
    write (10, *) '286 g440'
    goto 440
end if
!      IF (IBRNCH == 1 .AND. INFO /= 3) GOTO 440
!
!     Set the recently calculated function values in a column of DATMAT. This
!     array has a column for each vertex of the current simplex, the entries of
!     each column being the values of the constraint functions (if any)
!     followed by the objective function and the greatest constraint violation
!     at the vertex.
!
do K = 1, MPP
    DATMAT(K, JDROP) = CON(K)
end do

if (NFVALS > NP) then
    write (10, *) '302 g130'
    goto 130
end if
! IF we do not go to 130 but continue to below, then NFVALS <= NP.
! Thus NFVALS may be NP = N+1 > N.
!
!     Exchange the new vertex of the initial simplex with the optimal vertex if
!     necessary. Then, if the initial simplex is not complete, pick its next
!     vertex and calculate the function values there.
!
if (JDROP <= N) then
    if (DATMAT(MP, NP) <= F) then
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
        if (NFVALS <= N) X(JDROP) = SIM(JDROP, NP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
        SIM(JDROP, NP) = X(JDROP)
        do K = 1, MPP
            DATMAT(K, JDROP) = DATMAT(K, NP)
            DATMAT(K, NP) = CON(K)
        end do
        do K = 1, JDROP
            SIM(JDROP, K) = -RHO
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!              TEMP=0.0
            TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do I = K, JDROP
                TEMP = TEMP - SIMI(I, K)
            end do
            SIMI(JDROP, K) = TEMP
        end do
    end if
end if

! 120
if (NFVALS <= N) then
    JDROP = NFVALS
    X(JDROP) = X(JDROP) + RHO
    write (10, *) '367 g40'
    goto 40
end if
130 IBRNCH = 1
!
!     Identify the optimal vertex of the current simplex.
!
140 PHIMIN = DATMAT(MP, NP) + PARMU * DATMAT(MPP, NP)
NBEST = NP
do J = 1, N
    TEMP = DATMAT(MP, J) + PARMU * DATMAT(MPP, J)
    if (TEMP < PHIMIN) then
        NBEST = J
        PHIMIN = TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ELSE IF (TEMP .EQ. PHIMIN .AND. PARMU .EQ. 0.0) THEN
    else if (abs(TEMP - PHIMIN) <= ZERO .and. PARMU <= ZERO) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (DATMAT(MPP, J) < DATMAT(MPP, NBEST)) NBEST = J
    end if
end do
!
!     Switch the best vertex into pole position if it is not there already,
!     and also update SIM, SIMI and DATMAT.
!
if (NBEST <= N) then
    do I = 1, MPP
        TEMP = DATMAT(I, NP)
        DATMAT(I, NP) = DATMAT(I, NBEST)
        DATMAT(I, NBEST) = TEMP
    end do
    do I = 1, N
        TEMP = SIM(I, NBEST)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          SIM(I,NBEST)=0.0
        SIM(I, NBEST) = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SIM(I, NP) = SIM(I, NP) + TEMP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMPA=0.0
        TEMPA = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do K = 1, N
            SIM(I, K) = SIM(I, K) - TEMP
            TEMPA = TEMPA - SIMI(K, I)
        end do
        SIMI(NBEST, I) = TEMPA
    end do
end if

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 2021-05-30
if (INFO == 3) then
    write (10, *) '420 g600'
    goto 600
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     Make an error return if SIGI is a poor approximation to the inverse of
!     the leading N by N submatrix of SIG.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ERROR=0.0
ERROR = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do I = 1, N
    do J = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
!      IF (I .EQ. J) TEMP=TEMP-1.0
        TEMP = ZERO
        if (I == J) TEMP = TEMP - ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do K = 1, N
            TEMP = TEMP + SIMI(I, K) * SIM(K, J)
        end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  200 ERROR=AMAX1(ERROR,ABS(TEMP))
!      IF (ERROR .GT. 0.1) THEN
        ERROR = DMAX1(ERROR, DABS(TEMP))
    end do
end do
if (.not. (ERROR <= TENTH)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (IPRINT >= 1) print 210
210 format(/3X, 'Return from subroutine COBYLA because ', 'rounding errors are becoming damaging.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    INFO = 7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write (10, *) '457 g600'
    goto 600
end if
!
!     Calculate the coefficients of the linear approximations to the objective
!     and constraint functions, placing minus the objective function gradient
!     after the constraint gradients in the array A. The vector W is used for
!     working space.
!
! 220
do K = 1, MP
    CON(K) = -DATMAT(K, NP)
    do J = 1, N
        W(J) = DATMAT(K, J) + CON(K)
    end do
    do I = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      TEMP=0.0
        TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do J = 1, N
            TEMP = TEMP + W(J) * SIMI(J, I)
        end do
        if (K == MP) TEMP = -TEMP
        A(I, K) = TEMP
    end do
end do
!
!     Calculate the values of sigma and eta, and set IFLAG=0 if the current
!     simplex is not acceptable.
!
IFLAG = 1
PARSIG = ALPHA * RHO
PARETA = BETA * RHO
do J = 1, N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      WSIG=0.0
!      WETA=0.0
    WSIG = ZERO
    WETA = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do I = 1, N
        WSIG = WSIG + SIMI(J, I)**2
        WETA = WETA + SIM(I, J)**2
    end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      VSIG(J)=1.0/SQRT(WSIG)
!      VETA(J)=SQRT(WETA)
    VSIG(J) = ONE / DSQRT(WSIG)
    VETA(J) = DSQRT(WETA)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (VSIG(J) < PARSIG .or. VETA(J) > PARETA) IFLAG = 0
end do
!
!     If a new vertex is needed to improve acceptability, then decide which
!     vertex to drop from the simplex.
!
if (IBRNCH == 1 .or. IFLAG == 1) then
    write (10, *) '515 g370'
    goto 370
end if
JDROP = 0
TEMP = PARETA
do J = 1, N
    if (VETA(J) > TEMP) then
        JDROP = J
        TEMP = VETA(J)
    end if
end do
if (JDROP == 0) then
    do J = 1, N
        if (VSIG(J) < TEMP) then
            JDROP = J
            TEMP = VSIG(J)
        end if
    end do
end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190822: If VETA or VSIG become NaN due to rounding errors,
! JDROP may end up being 0. If we continue, then a Segmentation Fault
! will happen because we will read SIM(:, JDROP) and VSIG(JDROP).
if (JDROP == 0) then
    if (IPRINT >= 1) print 286
286 format(/3X, 'Return from subroutine COBYLA because ', 'rounding errors are becoming damaging.')
    INFO = 7
    write (10, *) '542 g600'
    goto 600
end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820: See the comments below line number 480
do I = 1, N
    XDROP(I) = SIM(I, NP) + SIM(I, JDROP) ! JDROP<NP is guaranteed
end do
do K = 1, MPP
    DATDROP(K) = DATMAT(K, JDROP)
end do
write (10, *) '530 bs'
call SAVEX(XDROP(1:N), DATDROP(1:MPP), XSAV(1:N, 1:NSMAX), DATSAV(1:MPP, 1:NSMAX), N, M, NSAV, NSMAX, CTOL)
write (10, *) '530 as'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Calculate the step to the new vertex and its sign.
!
TEMP = GAMMA * RHO * VSIG(JDROP)
do I = 1, N
    DX(I) = TEMP * SIMI(JDROP, I)
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      CVMAXP=0.0
!      CVMAXM=0.0
CVMAXP = ZERO
CVMAXM = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do K = 1, MP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SUM=0.0
    SUM = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do I = 1, N
        SUM = SUM + A(I, K) * DX(I)
    end do
    if (K < MP) then
        TEMP = DATMAT(K, NP)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          CVMAXP=AMAX1(CVMAXP,-SUM-TEMP)
!          CVMAXM=AMAX1(CVMAXM,SUM-TEMP)
        CVMAXP = DMAX1(CVMAXP, -SUM - TEMP)
        CVMAXM = DMAX1(CVMAXM, SUM - TEMP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end if
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      DXSIGN=1.0
!      IF (PARMU*(CVMAXP-CVMAXM) .GT. SUM+SUM) DXSIGN=-1.0
DXSIGN = ONE
if (PARMU * (CVMAXP - CVMAXM) > SUM + SUM) DXSIGN = -ONE
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
    SIM(I, JDROP) = DX(I)
    TEMP = TEMP + SIMI(JDROP, I) * DX(I)
end do
do I = 1, N
    SIMI(JDROP, I) = SIMI(JDROP, I) / TEMP
end do
do J = 1, N
    if (J /= JDROP) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=0.0
        TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do I = 1, N
            TEMP = TEMP + SIMI(J, I) * DX(I)
        end do
        do I = 1, N
            SIMI(J, I) = SIMI(J, I) - TEMP * SIMI(JDROP, I)
        end do
    end if
    X(J) = SIM(J, NP) + DX(J)
end do
write (10, *) '624 g40'
goto 40
!
!     Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO.
!
370 IZ = 1
IZDOTA = IZ + N * N
IVMC = IZDOTA + N
ISDIRN = IVMC + MP
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
376         format(/3X, 'Return from subroutine COBYLA because ', 'rounding errors are becoming damaging.')
            INFO = 7
            write (10, *) '648 g600'
            goto 600
        end if
    end do
end do
do J = 1, MP
    do I = 1, N
        if (is_nan(A(I, J))) then
            INFO = -3
            write (10, *) '657 g600'
            goto 600
        end if
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write (10, *) nfvals, 'bt'
!call TRSTLP(N, M, A, CON, RHO, DX, IFULL, IACT)!, W(IZ), W(IZDOTA), W(IVMC), W(ISDIRN), W(IDXNEW))!, W(IVMD:))
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!call TRSTLP(N, M, A, CON, RHO, DX, IFULL, IACT, W(IVMD:))
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

call TRSTLP(N, M, A, CON, RHO, DX, IFULL, IACT)
write (10, *) nfvals, 'at'
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
CON(MP) = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do K = 1, MP
    SUM = CON(K)
    do I = 1, N
        SUM = SUM - A(I, K) * DX(I)
    end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (K .LT. MP) RESNEW=AMAX1(RESNEW,SUM)
    if (K < MP) RESNEW = DMAX1(RESNEW, SUM)
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
PREREC = DATMAT(MPP, NP) - RESNEW
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (PREREC .GT. 0.0) BARMU=SUM/PREREC
!      IF (PARMU .LT. 1.5*BARMU) THEN
!          PARMU=2.0*BARMU
if (PREREC > ZERO) BARMU = SUM / PREREC
if (PARMU < 1.5D0 * BARMU) then
    PARMU = TWO * BARMU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (IPRINT >= 2) print 410, PARMU
410 format(/3X, 'Increase in PARMU to', 1PE13.6)
    PHI = DATMAT(MP, NP) + PARMU * DATMAT(MPP, NP)
    do J = 1, N
        TEMP = DATMAT(MP, J) + PARMU * DATMAT(MPP, J)
        if (TEMP < PHI) then
            write (10, *) '706 g140'
            goto 140
        end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          IF (TEMP .EQ. PHI .AND. PARMU .EQ. 0.0) THEN
        if (abs(TEMP - PHI) <= ZERO .and. PARMU <= ZERO) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (DATMAT(MPP, J) < DATMAT(MPP, NP)) then
                write (10, *) '714 g140'
                goto 140
            end if
        end if
    end do
end if
PREREM = PARMU * PREREC - SUM
!
!     Calculate the constraint and objective functions at x(*). Then find the
!     actual reduction in the merit function.
!
do I = 1, N
    X(I) = SIM(I, NP) + DX(I)
end do
IBRNCH = 1
write (10, *) '729 g40'
goto 40
440 VMOLD = DATMAT(MP, NP) + PARMU * DATMAT(MPP, NP)
VMNEW = F + PARMU * RESMAX
TRURED = VMOLD - VMNEW
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (PARMU .EQ. 0.0 .AND. F .EQ. DATMAT(MP,NP)) THEN
if (PARMU <= ZERO .and. abs(F - DATMAT(MP, NP)) <= ZERO) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PREREM = PREREC
    TRURED = DATMAT(MPP, NP) - RESMAX
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
RATIO = ZERO
if (TRURED <= ZERO) RATIO = ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
JDROP = 0
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
        JDROP = J
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
if (L > 0) JDROP = L
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
if (JDROP == 0) then
    do I = 1, N
        XDROP(I) = X(I)
    end do
    do K = 1, MPP
        DATDROP(K) = CONSAV(K)
    end do
else ! JDROP < NP is guaranteed
    do I = 1, N
        XDROP(I) = SIM(I, NP) + SIM(I, JDROP)
    end do
    do K = 1, MPP
        DATDROP(K) = DATMAT(K, JDROP)
    end do
end if
write (10, *) '820 bs'
call SAVEX(XDROP(1:N), DATDROP(1:MPP), XSAV(1:N, 1:NSMAX), DATSAV(1:MPP, 1:NSMAX), N, M, NSAV, NSMAX, CTOL)
write (10, *) '820 as'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (JDROP == 0) then
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
    SIM(I, JDROP) = DX(I)
    TEMP = TEMP + SIMI(JDROP, I) * DX(I)
end do
do I = 1, N
    SIMI(JDROP, I) = SIMI(JDROP, I) / TEMP
end do
do J = 1, N
    if (J /= JDROP) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          TEMP=0.0
        TEMP = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do I = 1, N
            TEMP = TEMP + SIMI(J, I) * DX(I)
        end do
        do I = 1, N
            SIMI(J, I) = SIMI(J, I) - TEMP * SIMI(JDROP, I)
        end do
    end if
end do
do K = 1, MPP
    DATMAT(K, JDROP) = CON(K)
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
550 if (IFLAG == 0) then
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
        do K = 1, MP
            CMIN = DATMAT(K, NP)
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
        if (abs(DENOM) <= ZERO) then  ! DENOM <= ZERO???  Is it nonnegative?
            PARMU = ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else if (CMAX - CMIN < PARMU * DENOM) then
            PARMU = (CMAX - CMIN) / DENOM
        end if
    end if
    if (IPRINT >= 2) print 580, RHO, PARMU
580 format(/3X, 'Reduction in RHO to', 1PE13.6, '  and PARMU =', 1PE13.6)
    if (IPRINT == 2) then
        print 70, NFVALS, DATMAT(MP, NP), DATMAT(MPP, NP), (SIM(I, NP), I=1, IPTEM)
        if (IPTEM < N) print 80, (X(I), I=IPTEMP, N)
    end if
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
590 format(/3X, 'Normal return from subroutine COBYLA')
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
600 do K = 1, M
    write (10, *) nfvals, 'be'
    CON(K) = CONSAV(K)
end do
PARMU = max(PARMU, 1.0D2)
if (NFVALS >= 2) then ! See the comments above for why NFVALS>2
    call ISBETTER(F, RESMAX, DATMAT(MP, NP), DATMAT(MPP, NP), PARMU, CTOL, BETTER)
    if (BETTER) then
        do I = 1, N
            X(I) = SIM(I, NP)
        end do
        F = DATMAT(MP, NP)
        RESMAX = DATMAT(MPP, NP)
        do K = 1, M
            CON(K) = DATMAT(K, NP)
        end do
    end if
    RESREF = RESMAX
    if (is_nan(RESREF)) RESREF = HUGENUM
    do J = 1, min(NP - 1, NFVALS - 2)
! See the comments above for why to check these J
        if (DATMAT(MPP, J) <= RESREF) then
            call ISBETTER(F, RESMAX, DATMAT(MP, J), DATMAT(MPP, J), PARMU, CTOL, BETTER)
            if (BETTER) then
                do I = 1, N
                    X(I) = SIM(I, J) + SIM(I, NP)
                end do
                F = DATMAT(MP, J)
                RESMAX = DATMAT(MPP, J)
                do K = 1, M
                    CON(K) = DATMAT(K, J)
                end do
            end if
        end if
    end do
end if
if (NSAV >= 1) then ! Do the following only if NSAV >= 1.
!          DO J = 1, NSAV
    do J = NSAV, 1, -1  ! We start with the most recent point
        if (DATSAV(MPP, J) <= RESREF) then
            call ISBETTER(F, RESMAX, DATSAV(MP, J), DATSAV(MPP, J), PARMU, CTOL, BETTER)
            if (BETTER) then
                do I = 1, N
                    X(I) = XSAV(I, J)
                end do
                F = DATSAV(MP, J)
                RESMAX = DATSAV(MPP, J)
                do K = 1, M
                    CON(K) = DATSAV(K, J)
                end do
            end if
        end if
    end do
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
620 if (IPRINT >= 1) then
    print 70, NFVALS, F, RESMAX, (X(I), I=1, IPTEM)
    if (IPTEM < N) print 80, (X(I), I=IPTEMP, N)
end if

close (10)
return
end subroutine cobylb

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820: See the comments below line number 480
subroutine SAVEX(XDROP, DATDROP, XSAV, DATSAV, N, M, NSAV, NSMAX, CTOL)
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
! subroutine with PARMU = -ONE. Indeed, PARMU can be any negative
! number. This is because, due to the implementation of ISBETTER,
! X dominates Y (i.e., X is better than Y with PARMU < 0)
! ==> X is better than Y with any PARMU >= 0,
! ==> X is better than Y regardless of PARMU.
! Fot this reason, it is sufficient to save all the "dropped" X that are not
! dominated by any vectors in XSAV (as we do in this subroutine),
! because the other X are always worse than some vector in XSAV.
!

! Generic modules
use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TENTH, HUGENUM, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
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
real(RP), intent(IN) :: datdrop(:)  ! m+2
real(RP), intent(IN) :: xdrop(:)  ! n

! In-outputs
integer(IK), intent(INOUT) :: Nsav
real(RP), intent(INOUT) :: Datsav(:, :)  ! (M+2, NSMAX)
real(RP), intent(INOUT) :: Xsav(:, :) ! (N, NSMAX)

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

if (NSMAX <= 0) return ! Do nothing if NSMAX=0

MP = M + 1
MPP = M + 2
PARMU = -ONE ! See the comments above for why PARMU = -1

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
    call ISBETTER(DATDROP(MP), DATDROP(MPP), DATSAV(MP, I), DATSAV(MPP, I), PARMU, CTOL, BETTER)
    if (BETTER) return
! If XDROP dominates XSAV(:, I), then increment NREMOVE by 1 and save
! I as IREMOVE(NREMOVE).
    call ISBETTER(DATSAV(MP, I), DATSAV(MPP, I), DATDROP(MP), DATDROP(MPP), PARMU, CTOL, BETTER)
    if (BETTER) then
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
J = 1 ! Index of IREMOVE
K = 1 ! Index of the new XSAV
do I = 1, NSAV ! Index of the old XSAV
    if (I == IREMOVE(J)) then
        J = J + 1 ! Do nothing but incrementing J by 1
    else ! Set XSAV(:, K) = XSAV(:, I)
        do L = 1, N
            XSAV(L, K) = XSAV(L, I)
        end do
        do L = 1, MPP
            DATSAV(L, K) = DATSAV(L, I)
        end do
        K = K + 1 ! Increment K by 1
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
    do L = 1, MPP
        DATSAV(L, NSAV) = DATDROP(L)
    end do
end if

end subroutine SAVEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun 20190820:
subroutine ISBETTER(F0, R0, F, R, PARMU, CTOL, BETTER)
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
! Generic modules
use consts_mod, only : RP, IK, ZERO, HALF, TENTH, HUGENUM, DEBUGGING, SRNLEN
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS, NAN_X, NAN_INF_F
use infnan_mod, only : is_nan, is_posinf
use debug_mod, only : errstop
use output_mod, only : retmssg, rhomssg, fmssg
use lina_mod, only : calquad, inprod
implicit none

real(RP), intent(IN) :: F0
real(RP), intent(IN) :: R0
real(RP), intent(IN) :: F
real(RP), intent(IN) :: R
real(RP), intent(IN) :: Parmu
real(RP), intent(IN) :: Ctol
logical, intent(OUT) :: Better

! Local variables
logical :: f0infnan
logical :: finfnan
logical :: fle
logical :: flt
logical :: r0infnan
logical :: rinfnan
logical :: rle
logical :: rlt

Better = .false.

! As values of F0, R0, F, and R, we regard Inf and NaN being equivalent values (they are equally bad).
f0infnan = is_nan(F0) .or. (F0 > HUGENUM)     ! F0 = Inf or NaN?
r0infnan = is_nan(R0) .or. (R0 > HUGENUM)     ! R0 = Inf or NaN?

BETTER = .false.

! As values of F0, R0, F, and R, we regard Inf and NaN being equivalent values (they are equally bad).
f0infnan = is_nan(F0) .or. (F0 > HUGENUM)     ! F0 = Inf or NaN?
r0infnan = is_nan(R0) .or. (R0 > HUGENUM)     ! R0 = Inf or NaN?
finfnan = is_nan(F) .or. (F > HUGENUM)     ! F = Inf or NaN?
rinfnan = is_nan(R) .or. (R > HUGENUM)     ! R  = Inf or NaN?

! When PARMU >= 0 and F + PARMU*R < F0 + PARMU*R0 and R < CTOL (feasible),
! then (F, R) is better than (F0, R0).
! Note that we should not set BETTER=FALSE even if this inequlity does not
! hold, because one or both of the two sides may be NaN.
if (PARMU >= ZERO .and. F + PARMU * R < F0 + PARMU * R0 .and. R < CTOL) then
    BETTER = .true.
end if

! If R < CTOL and F is not Inf or NaN while (R0 < CTOL) is false (may
! be because R0 is NaN), then (F, R) is better than (F0, R0). We prefer
! feasible points (i.e., constraint violation is less than CTOL) to
! insfeasible ones.
if (R < CTOL .and. .not. (R0 < CTOL) .and. .not. FINFNAN) then
    BETTER = .true.
end if

! If F0 or R0 is Inf/NaN while neither F nor R is Inf/NaN, then (F, R)
! is better than (F0, R0).
if ((F0INFNAN .or. R0INFNAN) .and. .not. (FINFNAN .or. RINFNAN)) then
    BETTER = .true.
end if

FLT = (F0INFNAN .and. (.not. FINFNAN)) .or. (F < F0) ! F < F0?
FLE = (F0INFNAN .and. FINFNAN) .or. (F <= F0) .or. FLT! F <= F0?
RLT = (R0INFNAN .and. (.not. RINFNAN)) .or. (R < R0) ! R < R0?
RLE = (R0INFNAN .and. RINFNAN) .or. (R <= R0) .or. RLT! R <= R0?

! If (F < F0 and R <= R0) or (F <= F0 and R < R0) in the sense defined
! above, the (F, R) is better than (F0, R0).
if ((FLT .and. RLE) .or. (FLE .and. RLT)) BETTER = .true.

! If one of F and R is -Inf and the other is not Inf/Nan, while neither
! F0 nor R0 is -Inf, then the (F, R) is better than (F0, R0).
if ((F < -HUGENUM) .and. (.not. RINFNAN) .and. (F0 >= -HUGENUM) .and. (R0 >= -HUGENUM)) BETTER = .true.
if ((R < -HUGENUM) .and. (.not. FINFNAN) .and. (F0 >= -HUGENUM) .and. (R0 >= -HUGENUM)) BETTER = .true.
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module cobylb_mod
