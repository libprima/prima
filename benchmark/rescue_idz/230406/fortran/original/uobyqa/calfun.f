      SUBROUTINE CALFUN (N,X,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),Y(10,10)
      DO 10 J=1,N
      Y(1,J)=1.0D0
   10 Y(2,J)=2.0D0*X(J)-1.0D0
      DO 20 I=2,N
      DO 20 J=1,N
   20 Y(I+1,J)=2.0D0*Y(2,J)*Y(I,J)-Y(I-1,J)
      F=0.0D0
      NP=N+1
      IW=1
      DO 40 I=1,NP
      SUM=0.0D0
      DO 30 J=1,N
   30 SUM=SUM+Y(I,J)
      SUM=SUM/DFLOAT(N)
      IF (IW .GT. 0) SUM=SUM+1.0/DFLOAT(I*I-2*I)
      IW=-IW
   40 F=F+SUM*SUM
      RETURN
      END
