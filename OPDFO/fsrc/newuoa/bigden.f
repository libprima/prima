      SUBROUTINE BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1  KNEW,D,W,VLAG,BETA,S,WVEC,PROD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!-----------------------!!!!!!
      USE DIRTY_TEMPORARY_MOD4POWELL_MOD!
      !!!!!!-----------------------!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION XOPT(*),XPT(NPT,N),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  W(*),VLAG(*),S(*),WVEC(NDIM,*),PROD(NDIM,*)
      DIMENSION DEN(9),DENEX(9),PAR(9), WSAVE(NPT), DOLD(N), vtmp(n)
C
C     N is the number of variables.
C     NPT is the number of interpolation equations.
C     XOPT is the best interpolation point so far.
C     XPT contains the coordinates of the current interpolation points.
C     BMAT provides the last N columns of H.
C     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     KOPT is the index of the optimal interpolation point.
C     KNEW is the index of the interpolation point that is going to be moved.
C     D will be set to the step from XOPT to the new point, and on entry it
C       should be the D that was calculated by the last call of BIGLAG. The
C       length of the initial D provides a trust region bound on the final D.
C     W will be set to Wcheck for the final choice of D.
C     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
C     BETA will be set to the value that will occur in the updating formula
C       when the KNEW-th interpolation point is moved to its new position.
C     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
C       for working space.
C
C     D is calculated in a way that should provide a denominator with a large
C     modulus in the updating formula when the KNEW-th interpolation point is
C     shifted to the new position XOPT+D.
C
C     Set some constants.
C
      !HALF=0.5D0
      !ONE=1.0D0
      !QUART=0.25D0
      !TWO=2.0D0
      !ZERO=0.0D0
      !TWOPI=8.0D0*DATAN(ONE)
      NPTM=NPT-N-1
C
C     Store the first NPT elements of the KNEW-th column of H in W(N+1)
C     to W(N+NPT).
C
      DO K=1,NPT
          W(N+K)=ZERO
      END DO
      DO J=1,NPTM
          TEMP=ZMAT(KNEW,J)
          IF (J < IDZ) TEMP=-TEMP
          DO K=1,NPT
              W(N+K)=W(N+K)+TEMP*ZMAT(K,J)
          END DO
      END DO
      ALPHA=W(N+KNEW)
C
C     The initial search direction D is taken from the last call of BIGLAG,
C     and the initial S is set below, usually to the direction from X_OPT
C     to X_KNEW, but a different direction to an interpolation point may
C     be chosen, in order to prevent S from being nearly parallel to D.
C
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DELTA = norm(D(1:N))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DD=ZERO
      DS=ZERO
      SS=ZERO
      XOPTSQ=ZERO
      DO I=1,N
          DD=DD+D(I)**2
          S(I)=XPT(KNEW,I)-XOPT(I)
          DS=DS+D(I)*S(I)
          SS=SS+S(I)**2
          XOPTSQ=XOPTSQ+XOPT(I)**2
      END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: With the original code, if DS, DD, or SS is
C NaN, KSAV will not get a value. This may cause Segmentation Fault
C because XPT(KSAV, :) will later be accessed.
C      IF (DS*DS .GT. 0.99D0*DD*SS) THEN
      !IF (.NOT. (DS*DS <= 0.99D0*DD*SS)) THEN
      IF (.NOT. (DS**2 <= 0.99D0*DD*SS)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          KSAV=KNEW
          !DTEST=DS*DS/SS
          DTEST=DS**2/SS
          DO K=1,NPT
              IF (K /= KOPT) THEN
                  DSTEMP=ZERO
                  SSTEMP=ZERO
                  DO I=1,N
                      DIFF=XPT(K,I)-XOPT(I)
                      DSTEMP=DSTEMP+D(I)*DIFF
                      SSTEMP=SSTEMP+DIFF*DIFF
                  END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun 2019-08-29: See the comments below line number 30
C              IF (DSTEMP*DSTEMP/SSTEMP .LT. DTEST) THEN
!                  IF (.NOT. (DSTEMP*DSTEMP/SSTEMP >= DTEST)) THEN
                  IF (.NOT. (DSTEMP**2/SSTEMP >= DTEST)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      KSAV=K
!                      DTEST=DSTEMP*DSTEMP/SSTEMP
                      DTEST=DSTEMP**2/SSTEMP
                      DS=DSTEMP
                      SS=SSTEMP
                  END IF
              END IF
          END DO
          DO I=1,N
              S(I)=XPT(KSAV,I)-XOPT(I)
          END DO
      END IF
      !SSDEN=DD*SS-DS*DS
      SSDEN=DD*SS-DS**2
      ITERC=0
      DENSAV=ZERO
C
C     Begin the iteration by overwriting S with a vector that has the
C     required length and direction.
C
   70 ITERC=ITERC+1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SSDEN=DD*SS-DS**2
      !IF (SSDEN <= max(SQRT(epsilon(0.0D0)), 1.0D-8)*DD*SS) GOTO 340
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      TEMP=ONE/DSQRT(SSDEN)
      XOPTD=ZERO
      XOPTS=ZERO
      !DO I=1,N
      !    S(I)=TEMP*(DD*S(I)-DS*D(I))
      !    XOPTD=XOPTD+XOPT(I)*D(I)
      !    XOPTS=XOPTS+XOPT(I)*S(I)
      !END DO
!      S(1:N)=S(1:N)-dot_product(S(1:n),
!     1 D(1:n)/norm(D(1:n)))*(D(1:N)/norm(D(1:N)))
      s(1:n) = s(1:n) - project(s(1:n), d(1:n))
      TLRNC = min(1.0D-1, max(epsilon(0.0D0)**(0.25D0), 1.0D-4))
      if (norm(S(1:n)) <= TLRNC*sqrt(SS)) then
          goto 340
      end if
      S(1:N) = (S(1:N)/norm(S(1:N)))*norm(D(1:N))
      if(abs(dot_product(d(1:n), s(1:n))) >=
     1 0.1D0*norm(d(1:n))*norm(s(1:n)) .or. norm(s(1:n)) >=2.0D0*DELTA)
     1 then
         GOTO  340
      end if
      DO I=1,N
          XOPTD=XOPTD+XOPT(I)*D(I)
          XOPTS=XOPTS+XOPT(I)*S(I)
      END DO
C
C     Set the coefficients of the first two terms of BETA.
C
      TEMPA=HALF*XOPTD*XOPTD
      TEMPB=HALF*XOPTS*XOPTS
      DEN(1)=DD*(XOPTSQ+HALF*DD)+TEMPA+TEMPB
      DEN(2)=TWO*XOPTD*DD
      DEN(3)=TWO*XOPTS*DD
      DEN(4)=TEMPA-TEMPB
      DEN(5)=XOPTD*XOPTS
      DO I=6,9
          DEN(I)=ZERO
      END DO
C
C     Put the coefficients of Wcheck in WVEC.
C
      DO K=1,NPT
          TEMPA=ZERO
          TEMPB=ZERO
          TEMPC=ZERO
          DO I=1,N
              TEMPA=TEMPA+XPT(K,I)*D(I)
              TEMPB=TEMPB+XPT(K,I)*S(I)
              TEMPC=TEMPC+XPT(K,I)*XOPT(I)
          END DO
!          WVEC(K,1)=QUART*(TEMPA*TEMPA+TEMPB*TEMPB)
          WVEC(K,1)=QUART*(TEMPA**2+TEMPB**2)
          WVEC(K,2)=TEMPA*TEMPC
          WVEC(K,3)=TEMPB*TEMPC
!          WVEC(K,4)=QUART*(TEMPA*TEMPA-TEMPB*TEMPB)
          WVEC(K,4)=QUART*(TEMPA**2-TEMPB**2)
          WVEC(K,5)=HALF*TEMPA*TEMPB
      END DO
      DO I=1,N
          IP=I+NPT
          WVEC(IP,1)=ZERO
          WVEC(IP,2)=D(I)
          WVEC(IP,3)=S(I)
          WVEC(IP,4)=ZERO
          WVEC(IP,5)=ZERO
      END DO
C
C     Put the coefficents of THETA*Wcheck in PROD.
C
      DO JC=1,5
          NW=NPT
          IF (JC == 2 .OR. JC == 3) NW=NDIM
          DO K=1,NPT
              PROD(K,JC)=ZERO
          END DO
          DO J=1,NPTM
              SUMM=ZERO
              DO K=1,NPT
                  SUMM=SUMM+ZMAT(K,J)*WVEC(K,JC)
              END DO
              IF (J < IDZ) SUMM=-SUMM
              DO K=1,NPT
                  PROD(K,JC)=PROD(K,JC)+SUMM*ZMAT(K,J)
              END DO
          END DO
          IF (NW == NDIM) THEN
              DO K=1,NPT
                  SUMM=ZERO
                  DO J=1,N
                      SUMM=SUMM+BMAT(K,J)*WVEC(NPT+J,JC)
                  END DO
                  PROD(K,JC)=PROD(K,JC)+SUMM
              END DO
          END IF
          DO J=1,N
              SUMM=ZERO
              DO I=1,NW
                  SUMM=SUMM+BMAT(I,J)*WVEC(I,JC)
              END DO
              PROD(NPT+J,JC)=SUMM
          END DO
      END DO
C
C     Include in DEN the part of BETA that depends on THETA.
C
      DO K=1,NDIM
          SUMM=ZERO
          DO I=1,5
              PAR(I)=HALF*PROD(K,I)*WVEC(K,I)
              SUMM=SUMM+PAR(I)
          END DO
          DEN(1)=DEN(1)-PAR(1)-SUMM
          TEMPA=PROD(K,1)*WVEC(K,2)+PROD(K,2)*WVEC(K,1)
          TEMPB=PROD(K,2)*WVEC(K,4)+PROD(K,4)*WVEC(K,2)
          TEMPC=PROD(K,3)*WVEC(K,5)+PROD(K,5)*WVEC(K,3)
          DEN(2)=DEN(2)-TEMPA-HALF*(TEMPB+TEMPC)
          DEN(6)=DEN(6)-HALF*(TEMPB-TEMPC)
          TEMPA=PROD(K,1)*WVEC(K,3)+PROD(K,3)*WVEC(K,1)
          TEMPB=PROD(K,2)*WVEC(K,5)+PROD(K,5)*WVEC(K,2)
          TEMPC=PROD(K,3)*WVEC(K,4)+PROD(K,4)*WVEC(K,3)
          DEN(3)=DEN(3)-TEMPA-HALF*(TEMPB-TEMPC)
          DEN(7)=DEN(7)-HALF*(TEMPB+TEMPC)
          TEMPA=PROD(K,1)*WVEC(K,4)+PROD(K,4)*WVEC(K,1)
          DEN(4)=DEN(4)-TEMPA-PAR(2)+PAR(3)
          TEMPA=PROD(K,1)*WVEC(K,5)+PROD(K,5)*WVEC(K,1)
          TEMPB=PROD(K,2)*WVEC(K,3)+PROD(K,3)*WVEC(K,2)
          DEN(5)=DEN(5)-TEMPA-HALF*TEMPB
          DEN(8)=DEN(8)-PAR(4)+PAR(5)
          TEMPA=PROD(K,4)*WVEC(K,5)+PROD(K,5)*WVEC(K,4)
          DEN(9)=DEN(9)-HALF*TEMPA
      END DO
C
C     Extend DEN so that it holds all the coefficients of DENOM.
C
      SUMM=ZERO
      DO I=1,5
          PAR(I)=HALF*PROD(KNEW,I)**2
          SUMM=SUMM+PAR(I)
      END DO
      DENEX(1)=ALPHA*DEN(1)+PAR(1)+SUMM
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,2)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,4)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,5)
      DENEX(2)=ALPHA*DEN(2)+TEMPA+TEMPB+TEMPC
      DENEX(6)=ALPHA*DEN(6)+TEMPB-TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,3)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,5)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,4)
      DENEX(3)=ALPHA*DEN(3)+TEMPA+TEMPB-TEMPC
      DENEX(7)=ALPHA*DEN(7)+TEMPB+TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,4)
      DENEX(4)=ALPHA*DEN(4)+TEMPA+PAR(2)-PAR(3)
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,5)
      DENEX(5)=ALPHA*DEN(5)+TEMPA+PROD(KNEW,2)*PROD(KNEW,3)
      DENEX(8)=ALPHA*DEN(8)+PAR(4)-PAR(5)
      DENEX(9)=ALPHA*DEN(9)+PROD(KNEW,4)*PROD(KNEW,5)
C
C     Seek the value of the angle that maximizes the modulus of DENOM.
C
      PAR(1)=ONE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !SUMM=DENEX(1)+DENEX(2)+DENEX(4)+DENEX(6)+DENEX(8)
      angle = 0.0D0 *((TWO*PI)/DFLOAT(50))
      par(2:8:2) = cos(angle * [1.0D0, 2.0D0, 3.0D0, 4.0D0])
      par(3:9:2) = sin(angle * [1.0D0, 2.0D0, 3.0D0, 4.0D0])
      SUMM = inprod(denex(1:9), par(1:9))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DENOLD=SUMM
      DENMAX=SUMM
      ISAVE=0
      IU=49
      TEMP=(TWO*PI)/DFLOAT(IU+1)
      DO I=1,IU
          ANGLE=DFLOAT(I)*TEMP
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !PAR(2)=DCOS(ANGLE)
          !PAR(3)=DSIN(ANGLE)
          !DO J=4,8,2
          !    PAR(J)=PAR(2)*PAR(J-2)-PAR(3)*PAR(J-1)
          !    PAR(J+1)=PAR(2)*PAR(J-1)+PAR(3)*PAR(J-2)
          !END DO
          par(2:8:2) = cos(angle * [1.0D0, 2.0D0, 3.0D0, 4.0D0])
          par(3:9:2) = sin(angle * [1.0D0, 2.0D0, 3.0D0, 4.0D0])
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SUMMOLD=SUMM
          SUMM=ZERO
          DO J=1,9
              SUMM=SUMM+DENEX(J)*PAR(J)
          END DO
          IF (DABS(SUMM) > DABS(DENMAX)) THEN
              DENMAX=SUMM
              ISAVE=I
              TEMPA=SUMMOLD
          ELSE IF (I == ISAVE+1) THEN
              TEMPB=SUMM
          END IF
      END DO
      IF (ISAVE == 0) TEMPA=SUMM
      IF (ISAVE == IU) TEMPB=DENOLD
      STEP=ZERO
      IF (TEMPA /= TEMPB) THEN
          TEMPA=TEMPA-DENMAX
          TEMPB=TEMPB-DENMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
C
C     Calculate the new parameters of the denominator, the new VLAG vector
C     and the new D. Then test for convergence.
C
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !PAR(2)=DCOS(ANGLE)
      !PAR(3)=DSIN(ANGLE)
      !DO J=4,8,2
      !    PAR(J)=PAR(2)*PAR(J-2)-PAR(3)*PAR(J-1)
      !    PAR(J+1)=PAR(2)*PAR(J-1)+PAR(3)*PAR(J-2)
      !END DO
      par(2:8:2) = cos(angle * [1.0D0, 2.0D0, 3.0D0, 4.0D0])
      par(3:9:2) = sin(angle * [1.0D0, 2.0D0, 3.0D0, 4.0D0])
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      BETA=ZERO
      DENMAX=ZERO
      DO J=1,9
          BETA=BETA+DEN(J)*PAR(J)
          DENMAX=DENMAX+DENEX(J)*PAR(J)
      END DO
      DO K=1,NDIM
          VLAG(K)=ZERO
          DO J=1,5
              VLAG(K)=VLAG(K)+PROD(K,J)*PAR(J)
          END DO
      END DO
      TAU=VLAG(KNEW)
      DD=ZERO
      TEMPA=ZERO
      TEMPB=ZERO
      dold(1:n) = d(1:n)
      DO I=1,N
          !D(I)=PAR(2)*D(I)+PAR(3)*S(I)
          D(I)=cos(angle)*D(I)+sin(angle)*S(I)
          W(I)=XOPT(I)+D(I)
          DD=DD+D(I)**2
          TEMPA=TEMPA+D(I)*W(I)
          TEMPB=TEMPB+W(I)*W(I)
      END DO
      ! Exit in case of Inf/NaN in D.
      temp = sum(abs(d(1:n)))
      if (temp /= temp .or. temp > huge(1.0D0)) then
         d(1:n) = dold(1:n)
         goto 340
      end if
      IF (ITERC >= N) GOTO 340
      IF (ITERC > 1) DENSAV=DMAX1(DENSAV,DENOLD)
      IF (DABS(DENMAX) <= 1.1D0*DABS(DENSAV)) GOTO 340
      DENSAV=DENMAX
C
C     Set S to half the gradient of the denominator with respect to D.
C     Then branch for the next iteration.
C
      DO I=1,N
          TEMP=TEMPA*XOPT(I)+TEMPB*D(I)-VLAG(NPT+I)
          S(I)=TAU*BMAT(KNEW,I)+ALPHA*TEMP
      END DO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      vtmp(1:n) = s(1:n) ; s(1:n) = zero
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO K=1,NPT
          SUMM=ZERO
          DO J=1,N
              SUMM=SUMM+XPT(K,J)*W(J)
          END DO
          TEMP=(TAU*W(N+K)-ALPHA*VLAG(K))*SUMM
          DO I=1,N
              S(I)=S(I)+TEMP*XPT(K,I)
          END DO
      END DO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      s(1:n) = vtmp(1:n) + s(1:n)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SS=ZERO
      DS=ZERO
      DO I=1,N
          SS=SS+S(I)**2
          DS=DS+D(I)*S(I)
      END DO
      !SSDEN=DD*SS-DS*DS
      !IF (SSDEN >= 1.0D-8*DD*SS) GOTO 70
      !SSDEN=DD*SS-DS**2
      !IF (SSDEN >= max(SQRT(epsilon(0.0D0)), 1.0D-8)*DD*SS) GOTO 70
      GOTO 70
C
C     Set the vector W before the RETURN from the subroutine.
C
  340 DO K=1,NDIM
          W(K)=ZERO
          DO J=1,5
              W(K)=W(K)+WVEC(K,J)*PAR(J)
          END DO
      END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      wsave = w(1:npt)
      w(1:npt) = matprod(xpt, d(1:n))
      w(1:npt)=w(1:npt)*(0.5D0*w(1:npt)+matprod(xpt,xopt(1:n)))
!      if (norm2(wsave-w(1:npt))/max(norm2(w(1:npt)), 1.0D0) > 1.0D-13)
!     &then
!        write(*,*) norm2(wsave-w(1:npt))/max(norm2(w(1:npt)), 1.0D0)
!        open(1,file ='w',status='old',position='append',action='write')
!        write(1,*) norm2(wsave-w(1:npt))/max(norm2(w(1:npt)), 1.0D0)
!        close(1)
!      end if
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      VLAG(KOPT)=VLAG(KOPT)+ONE
      RETURN
      END
