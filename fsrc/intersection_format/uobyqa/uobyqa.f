! This is the intersection-format version of uobyqa.f.
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


      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC&
     &CCCCC
      C SUBROUTINE UOBYQA (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
            SUBROUTINE UOBYQA (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,F,INFO,
           1 FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC&
     &CCCCC
      C IMPLICIT REAL*8 (A-H,O-Z)
            IMPLICIT REAL(KIND(0.0D0)) (A-H,O-Z)
            IMPLICIT INTEGER (I-N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DIMENSION X(*),W(*)
      C
      C This subroutine seeks the least value of a function of many vari&
     &ables,
      C by a trust region method that forms quadratic models by interpol&
     &ation.
      C The algorithm is described in "UOBYQA: unconstrained optimizatio&
     &n by
      C quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA1&
     &4,
      C University of Cambridge. The arguments of the subroutine are as &
     &follows.
      C
      C N must be set to the number of variables and must be at least tw&
     &o.
      C Initial values of the variables must be set in X(1),X(2),...,X(N&
     &). They
      C will be changed to the values that give the least calculated F.
      C RHOBEG and RHOEND must be set to the initial and final values of&
     &a trust
      C region radius, so both must be positive with RHOEND<=RHOBEG. Typ&
     &ically
      C RHOBEG should be about one tenth of the greatest expected change&
     &to a
      C variable, and RHOEND should indicate the accuracy that is requir&
     &ed in
      C the final values of the variables.
      C The value of IPRINT should be set to 0, 1, 2 or 3, which control&
     &s the
      C amount of printing. Specifically, there is no output if IPRINT=0&
     &and
      C there is output only at the return if IPRINT=1. Otherwise, each &
     &new
      C value of RHO is printed, with the best vector of variables so fa&
     &r and
      C the corresponding value of the objective function. Further, each&
     &new
      C value of F with its variables are output if IPRINT=3.
      C MAXFUN must be set to an upper bound on the number of calls of C&
     &ALFUN.
      C The array W will be used for working space. Its length must be a&
     &t least
      C ( N**4 + 8*N**3 + 23*N**2 + 42*N + max [ 2*N**2 + 4, 18*N ] ) / &
     &4.
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC&
     &CCCCC
      C F is the objective function value when the algorithm exit.
      C INFO is the exit flag, which can be set to:
      C 0: the lower bound for the trust region radius is reached.
      C 1: the target function value is reached.
      C 2: a trust region step has failed to reduce the quadratic model.
      C 3: the objective function has been evaluated MAXFUN times.
      C 4: much cancellation in a denominator.
      C 5: NPT is not in the required interval.
      C 6: one of the difference XU(I)-XL(I) is less than 2*RHOBEG.
      C 7: rounding errors are becoming damaging.
      C 8: rounding errors prevent reasonable changes to X.
      C 9: the denominator of the updating formule is zero.
      C 10: N should not be less than 2.
      C 11: MAXFUN is less than NPT+1.
      C 12: the gradient of constraint is zero.
      C -1: NaN occurs in x.
      C -2: the objective function returns a NaN or nearly infinite
      C value.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      C
      C SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must &
     &set F to
      C the value of the objective function for the variables X(1),X(2),&
     &...,X(N).
      C
      C Partition the working space array, so that different parts of it&
     &can be
      C treated separately by the subroutine that performs the main calc&
     &ulation.
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
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC&
     &CCCCC
      C Zaikun, 2020-05-05
      C When the data is passed from the interfaces to the Fortran code,&
     &RHOBEG,
      C and RHOEND may change a bit (due to rounding ???). It was oberve&
     &d in
      C a MATLAB test that MEX passed 1 to Fortran as 0.9999999999999997&
     &8.
      C If we set RHOEND = RHOBEG in the interfaces, then it may happen
      C that RHOEND > RHOBEG. That is why we do the following.
            RHOEND = MIN(RHOBEG, RHOEND)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CALL UOBYQB (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,NPT,W(IXB),W(IX&
     &O),
      CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC&
     &CCCCC
      C 1 W(IXN),W(IXP),W(IPQ),W(IPL),W(IH),W(IG),W(ID),W(IVL),W(IW))
           1 W(IXN),W(IXP),W(IPQ),W(IPL),W(IH),W(IG),W(ID),W(IVL),W(IW),&
     &F,
           2 INFO,FTARGET)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            RETURN
            END