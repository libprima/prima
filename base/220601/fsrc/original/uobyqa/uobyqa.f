      SUBROUTINE UOBYQA (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*)
C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     The algorithm is described in "UOBYQA: unconstrained optimization by
C     quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA14,
C     University of Cambridge. The arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
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
C       ( N**4 + 8*N**3 + 23*N**2 + 42*N + max [ 2*N**2 + 4, 18*N ] ) / 4.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NPT=(N*N+3*N+2)/2
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IPQ=IXP+N*NPT
      IPL=IPQ+NPT-1
      IH=IPL+(NPT-1)*NPT
      IG=IH+N*N
      ID=IG+N
      IVL=IH
      IW=ID+N
      CALL UOBYQB (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,NPT,W(IXB),W(IXO),
     1  W(IXN),W(IXP),W(IPQ),W(IPL),W(IH),W(IG),W(ID),W(IVL),W(IW))
      RETURN
      END
