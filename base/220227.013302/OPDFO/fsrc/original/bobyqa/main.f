C
C     Test problem for BOBYQA, the objective function being the sum of
C     the reciprocals of all pairwise distances between the points P_I,
C     I=1,2,...,M in two dimensions, where M=N/2 and where the components
C     of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables
C     defines the M points P_I. The initial X gives equally spaced points
C     on a circle. Four different choices of the pairs (N,NPT) are tried,
C     namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local
C     minimum that is not global occurs in both the N=10 cases. The details
C     of the results are highly sensitive to computer rounding errors. The
C     choice IPRINT=2 provides the current X and optimal F so far whenever
C     RHO is reduced. The bound constraints of the problem require every
C     component of X to be in the interval [-1,1].
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(100),XL(100),XU(100),W(500000)
      TWOPI=8.0D0*DATAN(1.0D0)
      BDL=-1.0D0
      BDU=1.0D0
      IPRINT=2
      MAXFUN=500000
      RHOBEG=1.0D-1
      RHOEND=1.0D-6
      M=5
   10 N=2*M
      DO 20 I=1,N
      XL(I)=BDL
   20 XU(I)=BDU
      DO 50 JCASE=1,2
      NPT=N+6
      IF (JCASE .EQ. 2) NPT=2*N+1
      PRINT 30, M,N,NPT
   30 FORMAT (//5X,'2D output with M =',I4,',  N =',I4,
     1  '  and  NPT =',I4)
      DO 40 J=1,M
      TEMP=DFLOAT(J)*TWOPI/DFLOAT(M)
      X(2*J-1)=DCOS(TEMP)
   40 X(2*J)=DSIN(TEMP)
      CALL BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
   50 CONTINUE
      M=M+M
      IF (M .LE. 10) GOTO 10
      STOP
      END
