C
C     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6,8.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(10),W(10000)
      IPRINT=2
      MAXFUN=5000
      RHOEND=1.0D-8
      DO 30 N=2,8,2
      DO 10 I=1,N
   10 X(I)=DFLOAT(I)/DFLOAT(N+1)
      RHOBEG=0.2D0*X(1)
      PRINT 20, N
   20 FORMAT (//5X,'******************'/5X,
     1  'Results with N =',I2,/5X,'******************')
      CALL UOBYQA (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
   30 CONTINUE
      STOP
      END
