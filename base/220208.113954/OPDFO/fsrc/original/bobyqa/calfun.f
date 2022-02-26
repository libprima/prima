      SUBROUTINE CALFUN (N,X,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*)
      F=0.0D0
      DO 10 I=4,N,2
      DO 10 J=2,I-2,2
      TEMP=(X(I-1)-X(J-1))**2+(X(I)-X(J))**2
      TEMP=DMAX1(TEMP,1.0D-6)
   10 F=F+1.0D0/DSQRT(TEMP)
      RETURN
      END
