      SUBROUTINE LAGMAX (N,G,H,RHO,D,V,VMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION G(*),H(N,*),D(*),V(*)
C
C     N is the number of variables of a quadratic objective function, Q say.
C     G is the gradient of Q at the origin.
C     H is the symmetric Hessian matrix of Q. Only the upper triangular and
C       diagonal parts need be set.
C     RHO is the trust region radius, and has to be positive.
C     D will be set to the calculated vector of variables.
C     The array V will be used for working space.
C     VMAX will be set to |Q(0)-Q(D)|.
C
C     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| .LEQ. RHO
C     requires of order N**3 operations, but sometimes it is adequate if
C     |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
C     subroutine provides such a solution in only of order N**2 operations,
C     where the claim of accuracy has been tested by numerical experiments.
C
C     Preliminary calculations.
C
      HALF=0.5D0
      HALFRT=DSQRT(HALF)
      ONE=1.0D0
      ZERO=0.0D0
C
C     Pick V such that ||HV|| / ||V|| is large.
C
      HMAX=ZERO
      DO 20 I=1,N
      SUM=ZERO
      DO 10 J=1,N
      H(J,I)=H(I,J)
   10 SUM=SUM+H(I,J)**2
      IF (SUM .GT. HMAX) THEN
          HMAX=SUM
          K=I
      END IF
   20 CONTINUE
      DO 30 J=1,N
   30 V(J)=H(K,J)
C
C     Set D to a vector in the subspace spanned by V and HV that maximizes
C     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
C     The vector that has the name D at label 60 used to be the vector W.
C
      VSQ=ZERO
      VHV=ZERO
      DSQ=ZERO
      DO 50 I=1,N
      VSQ=VSQ+V(I)**2
      D(I)=ZERO
      DO 40 J=1,N
   40 D(I)=D(I)+H(I,J)*V(J)
      VHV=VHV+V(I)*D(I)
   50 DSQ=DSQ+D(I)**2
      IF (VHV*VHV .LE. 0.9999D0*DSQ*VSQ) THEN
          TEMP=VHV/VSQ
          WSQ=ZERO
          DO 60 I=1,N
          D(I)=D(I)-TEMP*V(I)
   60     WSQ=WSQ+D(I)**2
          WHW=ZERO
          RATIO=DSQRT(WSQ/VSQ)
          DO 80 I=1,N
          TEMP=ZERO
          DO 70 J=1,N
   70     TEMP=TEMP+H(I,J)*D(J)
          WHW=WHW+TEMP*D(I)
   80     V(I)=RATIO*V(I)
          VHV=RATIO*RATIO*VHV
          VHW=RATIO*WSQ
          TEMP=HALF*(WHW-VHV)
          TEMP=TEMP+DSIGN(DSQRT(TEMP**2+VHW**2),WHW+VHV)
          DO 90 I=1,N
   90     D(I)=VHW*V(I)+TEMP*D(I)
      END IF
C
C     We now turn our attention to the subspace spanned by G and D. A multiple
C     of the current D is returned if that choice seems to be adequate.
C
      GG=ZERO
      GD=ZERO
      DD=ZERO
      DHD=ZERO
      DO 110 I=1,N
      GG=GG+G(I)**2
      GD=GD+G(I)*D(I)
      DD=DD+D(I)**2
      SUM=ZERO
      DO 100 J=1,N
  100 SUM=SUM+H(I,J)*D(J)
  110 DHD=DHD+SUM*D(I)
      TEMP=GD/GG
      VV=ZERO
      SCALE=DSIGN(RHO/DSQRT(DD),GD*DHD)
      DO 120 I=1,N
      V(I)=D(I)-TEMP*G(I)
      VV=VV+V(I)**2
  120 D(I)=SCALE*D(I)
      GNORM=DSQRT(GG)
      IF (GNORM*DD .LE. 0.5D-2*RHO*DABS(DHD) .OR.
     1  VV/DD .LE. 1.0D-4) THEN
          VMAX=DABS(SCALE*(GD+HALF*SCALE*DHD))
          GOTO 170
      END IF
C
C     G and V are now orthogonal in the subspace spanned by G and D. Hence
C     we generate an orthonormal basis of this subspace such that (D,HV) is
C     negligible or zero, where D and V will be the basis vectors.
C
      GHG=ZERO
      VHG=ZERO
      VHV=ZERO
      DO 140 I=1,N
      SUM=ZERO
      SUMV=ZERO
      DO 130 J=1,N
      SUM=SUM+H(I,J)*G(J)
  130 SUMV=SUMV+H(I,J)*V(J)
      GHG=GHG+SUM*G(I)
      VHG=VHG+SUMV*G(I)
  140 VHV=VHV+SUMV*V(I)
      VNORM=DSQRT(VV)
      GHG=GHG/GG
      VHG=VHG/(VNORM*GNORM)
      VHV=VHV/VV
      IF (DABS(VHG) .LE. 0.01D0*DMAX1(DABS(GHG),DABS(VHV))) THEN
          VMU=GHG-VHV
          WCOS=ONE
          WSIN=ZERO
      ELSE
          TEMP=HALF*(GHG-VHV)
          VMU=TEMP+DSIGN(DSQRT(TEMP**2+VHG**2),TEMP)
          TEMP=DSQRT(VMU**2+VHG**2)
          WCOS=VMU/TEMP
          WSIN=VHG/TEMP
      END IF
      TEMPA=WCOS/GNORM
      TEMPB=WSIN/VNORM
      TEMPC=WCOS/VNORM
      TEMPD=WSIN/GNORM
      DO 150 I=1,N
      D(I)=TEMPA*G(I)+TEMPB*V(I)
  150 V(I)=TEMPC*V(I)-TEMPD*G(I)
C
C     The final D is a multiple of the current D, V, D+V or D-V. We make the
C     choice from these possibilities that is optimal.
C
      DLIN=WCOS*GNORM/RHO
      VLIN=-WSIN*GNORM/RHO
      TEMPA=DABS(DLIN)+HALF*DABS(VMU+VHV)
      TEMPB=DABS(VLIN)+HALF*DABS(GHG-VMU)
      TEMPC=HALFRT*(DABS(DLIN)+DABS(VLIN))+0.25D0*DABS(GHG+VHV)
      IF (TEMPA .GE. TEMPB .AND. TEMPA .GE. TEMPC) THEN
          TEMPD=DSIGN(RHO,DLIN*(VMU+VHV))
          TEMPV=ZERO
      ELSE IF (TEMPB .GE. TEMPC) THEN
          TEMPD=ZERO
          TEMPV=DSIGN(RHO,VLIN*(GHG-VMU))
      ELSE
          TEMPD=DSIGN(HALFRT*RHO,DLIN*(GHG+VHV))
          TEMPV=DSIGN(HALFRT*RHO,VLIN*(GHG+VHV))
      END IF
      DO 160 I=1,N
  160 D(I)=TEMPD*D(I)+TEMPV*V(I)
      VMAX=RHO*RHO*DMAX1(TEMPA,TEMPB,TEMPC)
  170 RETURN
      END
