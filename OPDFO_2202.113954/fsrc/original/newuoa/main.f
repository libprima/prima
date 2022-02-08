C
C     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
C     with NPT = 2N+1.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(10),W(10000)
      IPRINT=2
      MAXFUN=5000
      RHOEND=1.0D-6
      DO 30 N=2,8,2
      NPT=2*N+1
      DO 10 I=1,N
   10 X(I)=DFLOAT(I)/DFLOAT(N+1)
      RHOBEG=0.2D0*X(1)
      PRINT 20, N,NPT
   20 FORMAT (//4X,'Results with N =',I2,' and NPT =',I3)
      CALL NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
   30 CONTINUE
      STOP
      END
