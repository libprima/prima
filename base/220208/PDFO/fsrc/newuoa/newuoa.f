CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,F,INFO,
     1  FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIMENSION X(*),W(*)
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     F is the objective function value when the algorithm exit.
C     INFO is the exit flag, which can be set to:
C       0: the lower bound for the trust region radius is reached.
C       1: the target function value is reached.
C       2: a trust region step has failed to reduce the quadratic model.
C       3: the objective function has been evaluated MAXFUN times.
C       4: much cancellation in a denominator.
C       5: NPT is not in the required interval.
C       6: one of the difference XU(I)-XL(I) is less than 2*RHOBEG.
C       7: rounding errors are becoming damaging.
C       8: rounding errors prevent reasonable changes to X.
C       9: the denominator of the updating formule is zero.
C       10: N should not be less than 2.
C       11: MAXFUN is less than NPT+1.
C       12: the gradient of constraint is zero.
C       -1: NaN occurs in x.
C       -2: the objective function returns a NaN or nearly infinite
C           value.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NP=N+1
      NPTM=NPT-NP
      IF (NPT < N+2 .OR. NPT > ((N+2)*NP)/2) THEN
          IF (IPRINT > 0) PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      ' the required interval')
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          INFO=5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IFV=IXP+N*NPT
      IGQ=IFV+NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
C
C     The above settings provide a partition of W for subroutine NEWUOB.
C     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
C     W plus the space that is needed by the last array of NEWUOB.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Zaikun, 2020-05-05
C When the data is passed from the interfaces to the Fortran code, RHOBEG, 
C and RHOEND may change a bit (due to rounding ???). It was oberved in
C a MATLAB test that MEX passed 1 to Fortran as 0.99999999999999978.
C If we set RHOEND = RHOBEG in the interfaces, then it may happen
C that RHOEND > RHOBEG. That is why we do the following. 
      RHOEND = MIN(RHOBEG, RHOEND)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW))
     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW),F,INFO,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   20 RETURN
      END
