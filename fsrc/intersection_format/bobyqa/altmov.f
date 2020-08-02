! This is the intersection-format version of altmov.f.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection format, each continued line has an ampersand at
! column 73, and each continuation line has an ampersand at column 6.
! A Fortran file in such a format can be compiled both as a fixed-format
! file and as a free-format file.
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 03-Aug-2020.


            SUBROUTINE ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
           1 KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,GLAG,HCOL,W)
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC&
     &CCCCC
      C IMPLICIT REAL*8 (A-H,O-Z)
            IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
            IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DIMENSION XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),
           1 SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
      C
      C The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all&
     &have
      C the same meanings as the corresponding arguments of BOBYQB.
      C KOPT is the index of the optimal interpolation point.
      C KNEW is the index of the interpolation point that is going to be&
     &moved.
      C ADELT is the current trust region bound.
      C XNEW will be set to a suitable new position for the interpolatio&
     &n point
      C XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust reg&
     &ion
      C bounds and it should provide a large denominator in the next cal&
     &l of
      C UPDATE. The step XNEW-XOPT from XOPT is restricted to moves alon&
     &g the
      C straight lines through XOPT and another interpolation point.
      C XALT also provides a large value of the modulus of the KNEW-th L&
     &agrange
      C function subject to the constraints that have been mentioned, it&
     &s main
      C difference from XNEW being that XALT-XOPT is a constrained versi&
     &on of
      C the Cauchy step within the trust region. An exception is that XA&
     &LT is
      C not calculated if all components of GLAG (see below) are zero.
      C ALPHA will be set to the KNEW-th diagonal element of the H matri&
     &x.
      C CAUCHY will be set to the square of the KNEW-th Lagrange functio&
     &n at
      C the step XALT-XOPT from XOPT for the vector XALT that is returne&
     &d,
      C except that CAUCHY is set to zero if XALT is not calculated.
      C GLAG is a working space vector of length N for the gradient of t&
     &he
      C KNEW-th Lagrange function at XOPT.
      C HCOL is a working space vector of length NPT for the second deri&
     &vative
      C coefficients of the KNEW-th Lagrange function.
      C W is a working space vector of length 2N that is going to hold t&
     &he
      C constrained Cauchy step from XOPT of the Lagrange function, foll&
     &owed
      C by the downhill version of XALT when the uphill step is calculat&
     &ed.
      C
      C Set the first NPT components of W to the leading elements of the
      C KNEW-th column of the H matrix.
      C
            HALF=0.5D0
            ONE=1.0D0
            ZERO=0.0D0
            CONST=ONE+DSQRT(2.0D0)
            DO K=1,NPT
                HCOL(K)=ZERO
            END DO
            DO J=1,NPT-N-1
                TEMP=ZMAT(KNEW,J)
                DO K=1,NPT
                    HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
                END DO
            END DO
            ALPHA=HCOL(KNEW)
            HA=HALF*ALPHA
      C
      C Calculate the gradient of the KNEW-th Lagrange function at XOPT.
      C
            DO I=1,N
                GLAG(I)=BMAT(KNEW,I)
            END DO
            DO K=1,NPT
                TEMP=ZERO
                DO J=1,N
                    TEMP=TEMP+XPT(K,J)*XOPT(J)
                END DO
                TEMP=HCOL(K)*TEMP
                DO I=1,N
                    GLAG(I)=GLAG(I)+TEMP*XPT(K,I)
                END DO
            END DO
      C
      C Search for a large denominator along the straight lines through &
     &XOPT
      C and another interpolation point. SLBD and SUBD will be lower and&
     &upper
      C bounds on the step along each of these lines in turn. PREDSQ wil&
     &l be
      C set to the square of the predicted denominator for each line. PR&
     &ESAV
      C will be set to the largest admissible value of PREDSQ that occur&
     &s.
      C
            PRESAV=ZERO
            DO K=1,NPT
                IF (K == KOPT) CYCLE
                DDERIV=ZERO
                DISTSQ=ZERO
                DO I=1,N
                    TEMP=XPT(K,I)-XOPT(I)
                    DDERIV=DDERIV+GLAG(I)*TEMP
                    DISTSQ=DISTSQ+TEMP*TEMP
                END DO
                SUBD=ADELT/DSQRT(DISTSQ)
                SLBD=-SUBD
                ILBD=0
                IUBD=0
                SUMIN=DMIN1(ONE,SUBD)
      C
      C Revise SLBD and SUBD if necessary because of the bounds in SL an&
     &d SU.
      C
                DO I=1,N
                    TEMP=XPT(K,I)-XOPT(I)
                    IF (TEMP > ZERO) THEN
                        IF (SLBD*TEMP < SL(I)-XOPT(I)) THEN
                            SLBD=(SL(I)-XOPT(I))/TEMP
                            ILBD=-I
                        END IF
                        IF (SUBD*TEMP > SU(I)-XOPT(I)) THEN
                            SUBD=DMAX1(SUMIN,(SU(I)-XOPT(I))/TEMP)
                            IUBD=I
                        END IF
                    ELSE IF (TEMP < ZERO) THEN
                        IF (SLBD*TEMP > SU(I)-XOPT(I)) THEN
                            SLBD=(SU(I)-XOPT(I))/TEMP
                            ILBD=I
                        END IF
                        IF (SUBD*TEMP < SL(I)-XOPT(I)) THEN
                            SUBD=DMAX1(SUMIN,(SL(I)-XOPT(I))/TEMP)
                            IUBD=-I
                        END IF
                    END IF
                END DO
      C
      C Seek a large modulus of the KNEW-th Lagrange function when the i&
     &ndex
      C of the other interpolation point on the line through XOPT is KNE&
     &W.
      C
                IF (K == KNEW) THEN
                    DIFF=DDERIV-ONE
                    STEP=SLBD
                    VLAG=SLBD*(DDERIV-SLBD*DIFF)
                    ISBD=ILBD
                    TEMP=SUBD*(DDERIV-SUBD*DIFF)
                    IF (DABS(TEMP) > DABS(VLAG)) THEN
                        STEP=SUBD
                        VLAG=TEMP
                        ISBD=IUBD
                    END IF
                    TEMPD=HALF*DDERIV
                    TEMPA=TEMPD-DIFF*SLBD
                    TEMPB=TEMPD-DIFF*SUBD
                    IF (TEMPA*TEMPB < ZERO) THEN
                        TEMP=TEMPD*TEMPD/DIFF
                        IF (DABS(TEMP) > DABS(VLAG)) THEN
                            STEP=TEMPD/DIFF
                            VLAG=TEMP
                            ISBD=0
                        END IF
                    END IF
      C
      C Search along each of the other lines through XOPT and another po&
     &int.
      C
                ELSE
                    STEP=SLBD
                    VLAG=SLBD*(ONE-SLBD)
                    ISBD=ILBD
                    TEMP=SUBD*(ONE-SUBD)
                    IF (DABS(TEMP) > DABS(VLAG)) THEN
                        STEP=SUBD
                        VLAG=TEMP
                        ISBD=IUBD
                    END IF
                    IF (SUBD > HALF) THEN
                        IF (DABS(VLAG) < 0.25D0) THEN
                            STEP=HALF
                            VLAG=0.25D0
                            ISBD=0
                        END IF
                    END IF
                    VLAG=VLAG*DDERIV
                END IF
      C
      C Calculate PREDSQ for the current line search and maintain PRESAV&
     &.
      C
                TEMP=STEP*(ONE-STEP)*DISTSQ
                PREDSQ=VLAG*VLAG*(VLAG*VLAG+HA*TEMP*TEMP)
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC&
     &CCCCC
      C Zaikun 2019-08-29: With the original code, if either PREDSQ or P&
     &RESAV
      C is NaN, KSAV/STPSAV/IBDSAV will not get a value. This may cause
      C Segmentation Fault.
      C IF (PREDSQ .GT. PRESAV) THEN
                IF (.NOT. (PREDSQ <= PRESAV)) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    PRESAV=PREDSQ
                    KSAV=K
                    STPSAV=STEP
                    IBDSAV=ISBD
                END IF
            END DO
      C
      C Construct XNEW in a way that satisfies the bound constraints exa&
     &ctly.
      C
            DO I=1,N
                TEMP=XOPT(I)+STPSAV*(XPT(KSAV,I)-XOPT(I))
                XNEW(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
            END DO
            IF (IBDSAV < 0) XNEW(-IBDSAV)=SL(-IBDSAV)
            IF (IBDSAV > 0) XNEW(IBDSAV)=SU(IBDSAV)
      C
      C Prepare for the iterative method that assembles the constrained &
     &Cauchy
      C step in W. The sum of squares of the fixed components of W is fo&
     &rmed in
      C WFIXSQ, and the free components of W are set to BIGSTP.
      C
            BIGSTP=ADELT+ADELT
            IFLAG=0
        100 WFIXSQ=ZERO
            GGFREE=ZERO
            DO I=1,N
                W(I)=ZERO
                TEMPA=DMIN1(XOPT(I)-SL(I),GLAG(I))
                TEMPB=DMAX1(XOPT(I)-SU(I),GLAG(I))
                IF (TEMPA > ZERO .OR. TEMPB < ZERO) THEN
                    W(I)=BIGSTP
                    GGFREE=GGFREE+GLAG(I)**2
                END IF
            END DO
            IF (GGFREE == ZERO) THEN
                CAUCHY=ZERO
                GOTO 200
            END IF
      C
      C Investigate whether more components of W can be fixed.
      C
        120 TEMP=ADELT*ADELT-WFIXSQ
            IF (TEMP > ZERO) THEN
                WSQSAV=WFIXSQ
                STEP=DSQRT(TEMP/GGFREE)
                GGFREE=ZERO
                DO I=1,N
                    IF (W(I) == BIGSTP) THEN
                        TEMP=XOPT(I)-STEP*GLAG(I)
                        IF (TEMP <= SL(I)) THEN
                            W(I)=SL(I)-XOPT(I)
                            WFIXSQ=WFIXSQ+W(I)**2
                        ELSE IF (TEMP >= SU(I)) THEN
                            W(I)=SU(I)-XOPT(I)
                            WFIXSQ=WFIXSQ+W(I)**2
                        ELSE
                            GGFREE=GGFREE+GLAG(I)**2
                        END IF
                    END IF
                END DO
                IF (WFIXSQ > WSQSAV .AND. GGFREE > ZERO) GOTO 120
            END IF
      C
      C Set the remaining free components of W and all components of XAL&
     &T,
      C except that W may be scaled later.
      C
            GW=ZERO
            DO I=1,N
                IF (W(I) == BIGSTP) THEN
                    W(I)=-STEP*GLAG(I)
                    XALT(I)=DMAX1(SL(I),DMIN1(SU(I),XOPT(I)+W(I)))
                ELSE IF (W(I) == ZERO) THEN
                    XALT(I)=XOPT(I)
                ELSE IF (GLAG(I) > ZERO) THEN
                    XALT(I)=SL(I)
                ELSE
                    XALT(I)=SU(I)
                END IF
                GW=GW+GLAG(I)*W(I)
            END DO
      C
      C Set CURV to the curvature of the KNEW-th Lagrange function along&
     &W.
      C Scale W by a factor less than one if that can reduce the modulus&
     &of
      C the Lagrange function at XOPT+W. Set CAUCHY to the final value o&
     &f
      C the square of this function.
      C
            CURV=ZERO
            DO K=1,NPT
                TEMP=ZERO
                DO J=1,N
                    TEMP=TEMP+XPT(K,J)*W(J)
                END DO
                CURV=CURV+HCOL(K)*TEMP*TEMP
            END DO
            IF (IFLAG == 1) CURV=-CURV
            IF (CURV > -GW .AND. CURV < -CONST*GW) THEN
                SCALE=-GW/CURV
                DO I=1,N
                    TEMP=XOPT(I)+SCALE*W(I)
                    XALT(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
                END DO
                CAUCHY=(HALF*GW*SCALE)**2
            ELSE
                CAUCHY=(GW+HALF*CURV)**2
            END IF
      C
      C If IFLAG is zero, then XALT is calculated as before after revers&
     &ing
      C the sign of GLAG. Thus two XALT vectors become available. The on&
     &e that
      C is chosen is the one that gives the larger value of CAUCHY.
      C
            IF (IFLAG == 0) THEN
                DO I=1,N
                    GLAG(I)=-GLAG(I)
                    W(N+I)=XALT(I)
                END DO
                CSAVE=CAUCHY
                IFLAG=1
                GOTO 100
            END IF
            IF (CSAVE > CAUCHY) THEN
                DO I=1,N
                    XALT(I)=W(N+I)
                END DO
                CAUCHY=CSAVE
            END IF
        200 RETURN
            END