      SUBROUTINE LAGMAX (N,G,H,RHO,D,V,VMAX)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use, non_intrinsic :: consts_mod, only : RP, IK
      implicit real(RP) (a - h, o - z)
      implicit integer(IK) (i - n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      HALFRT=SQRT(HALF)
      ONE=1.0D0
      ZERO=0.0D0
C
C     Pick V such that ||HV|| / ||V|| is large.
C
      HMAX=ZERO
      DO I=1,N
          SUM=ZERO
          DO J=1,N
              H(J,I)=H(I,J)
              SUM=SUM+H(I,J)**2
          END DO
          IF (SUM > HMAX) THEN
              HMAX=SUM
              K=I
          END IF
      END DO
      DO J=1,N
          V(J)=H(K,J)
      END DO
C
C     Set D to a vector in the subspace spanned by V and HV that maximizes
C     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
C     The vector that has the name D at label 60 used to be the vector W.
C
      VSQ=ZERO
      VHV=ZERO
      DSQ=ZERO
      DO I=1,N
          VSQ=VSQ+V(I)**2
          D(I)=ZERO
          DO J=1,N
              D(I)=D(I)+H(I,J)*V(J)
          END DO
          VHV=VHV+V(I)*D(I)
          DSQ=DSQ+D(I)**2
      END DO
      IF (VHV*VHV <= 0.9999D0*DSQ*VSQ) THEN
          TEMP=VHV/VSQ
          WSQ=ZERO
          DO I=1,N
              D(I)=D(I)-TEMP*V(I)
              WSQ=WSQ+D(I)**2
          END DO
          WHW=ZERO
          RATIO=SQRT(WSQ/VSQ)
          DO I=1,N
              TEMP=ZERO
              DO J=1,N
                  TEMP=TEMP+H(I,J)*D(J)
              END DO
              WHW=WHW+TEMP*D(I)
              V(I)=RATIO*V(I)
          END DO
          VHV=RATIO*RATIO*VHV
          VHW=RATIO*WSQ
          TEMP=HALF*(WHW-VHV)
          TEMP=TEMP+SIGN(SQRT(TEMP**2+VHW**2),WHW+VHV)
          DO I=1,N
              D(I)=VHW*V(I)+TEMP*D(I)
          END DO
      END IF
C
C     We now turn our attention to the subspace spanned by G and D. A multiple
C     of the current D is returned if that choice seems to be adequate.
C
      GG=ZERO
      GD=ZERO
      DD=ZERO
      DHD=ZERO
      DO I=1,N
          GG=GG+G(I)**2
          GD=GD+G(I)*D(I)
          DD=DD+D(I)**2
          SUM=ZERO
          DO J=1,N
              SUM=SUM+H(I,J)*D(J)
          END DO
          DHD=DHD+SUM*D(I)
      END DO
      TEMP=GD/GG
      VV=ZERO
      SCALE=SIGN(RHO/SQRT(DD),GD*DHD)
      DO I=1,N
          V(I)=D(I)-TEMP*G(I)
          VV=VV+V(I)**2
          D(I)=SCALE*D(I)
      END DO
      GNORM=SQRT(GG)
      IF (GNORM*DD <= 0.5D-2*RHO*ABS(DHD) .OR.
     1  VV/DD <= 1.0D-4) THEN
          VMAX=ABS(SCALE*(GD+HALF*SCALE*DHD))
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
      DO I=1,N
          SUM=ZERO
          SUMV=ZERO
          DO J=1,N
              SUM=SUM+H(I,J)*G(J)
              SUMV=SUMV+H(I,J)*V(J)
          END DO
          GHG=GHG+SUM*G(I)
          VHG=VHG+SUMV*G(I)
          VHV=VHV+SUMV*V(I)
      END DO
      VNORM=SQRT(VV)
      GHG=GHG/GG
      VHG=VHG/(VNORM*GNORM)
      VHV=VHV/VV
      IF (ABS(VHG) <= 0.01D0*MAX(ABS(GHG),ABS(VHV))) THEN
          VMU=GHG-VHV
          WCOS=ONE
          WSIN=ZERO
      ELSE
          TEMP=HALF*(GHG-VHV)
          VMU=TEMP+SIGN(SQRT(TEMP**2+VHG**2),TEMP)
          TEMP=SQRT(VMU**2+VHG**2)
          WCOS=VMU/TEMP
          WSIN=VHG/TEMP
      END IF
      TEMPA=WCOS/GNORM
      TEMPB=WSIN/VNORM
      TEMPC=WCOS/VNORM
      TEMPD=WSIN/GNORM
      DO I=1,N
          D(I)=TEMPA*G(I)+TEMPB*V(I)
          V(I)=TEMPC*V(I)-TEMPD*G(I)
      END DO
C
C     The final D is a multiple of the current D, V, D+V or D-V. We make the
C     choice from these possibilities that is optimal.
C
      DLIN=WCOS*GNORM/RHO
      VLIN=-WSIN*GNORM/RHO
      TEMPA=ABS(DLIN)+HALF*ABS(VMU+VHV)
      TEMPB=ABS(VLIN)+HALF*ABS(GHG-VMU)
      TEMPC=HALFRT*(ABS(DLIN)+ABS(VLIN))+0.25D0*ABS(GHG+VHV)
      IF (TEMPA >= TEMPB .AND. TEMPA >= TEMPC) THEN
          TEMPD=SIGN(RHO,DLIN*(VMU+VHV))
          TEMPV=ZERO
      ELSE IF (TEMPB >= TEMPC) THEN
          TEMPD=ZERO
          TEMPV=SIGN(RHO,VLIN*(GHG-VMU))
      ELSE
          TEMPD=SIGN(HALFRT*RHO,DLIN*(GHG+VHV))
          TEMPV=SIGN(HALFRT*RHO,VLIN*(GHG+VHV))
      END IF
      DO I=1,N
          D(I)=TEMPD*D(I)+TEMPV*V(I)
      END DO
      VMAX=RHO*RHO*MAX(TEMPA,TEMPB,TEMPC)
  170 RETURN
      END
