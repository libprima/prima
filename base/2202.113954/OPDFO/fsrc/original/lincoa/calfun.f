      SUBROUTINE CALFUN (N,X,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON FMAX
      DIMENSION X(*)
      ZERO=0.0D0
      F=FMAX
      V12=X(1)*X(5)-X(4)*X(2)
      V13=X(1)*X(8)-X(7)*X(2)
      V14=X(1)*X(11)-X(10)*X(2)
      V23=X(4)*X(8)-X(7)*X(5)
      V24=X(4)*X(11)-X(10)*X(5)
      V34=X(7)*X(11)-X(10)*X(8)
      DEL1=V23*X(12)-V24*X(9)+V34*X(6)
      IF (DEL1 .LE. ZERO) GOTO 10
      DEL2=-V34*X(3)-V13*X(12)+V14*X(9)
      IF (DEL2 .LE. ZERO) GOTO 10
      DEL3=-V14*X(6)+V24*X(3)+V12*X(12)
      IF (DEL3 .LE. ZERO) GOTO 10
      DEL4=-V12*X(9)+V13*X(6)-V23*X(3)
      IF (DEL4 .LE. ZERO) GOTO 10
      TEMP=(DEL1+DEL2+DEL3+DEL4)**3/(DEL1*DEL2*DEL3*DEL4)
      F=DMIN1(TEMP/6.0D0,FMAX)
   10 CONTINUE
      RETURN
      END
