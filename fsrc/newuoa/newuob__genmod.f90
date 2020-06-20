        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun 20 19:12:38 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NEWUOB__genmod
          INTERFACE 
            SUBROUTINE NEWUOB(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE,&
     &XOPT,XNEW,XPT,FVAL,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W,F,INFO,FTARGET&
     &)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(INOUT) :: X(N)
              REAL(KIND=8), INTENT(IN) :: RHOBEG
              REAL(KIND=8), INTENT(IN) :: RHOEND
              INTEGER(KIND=4), INTENT(IN) :: IPRINT
              INTEGER(KIND=4), INTENT(IN) :: MAXFUN
              REAL(KIND=8), INTENT(INOUT) :: XBASE(N)
              REAL(KIND=8), INTENT(INOUT) :: XOPT(N)
              REAL(KIND=8), INTENT(INOUT) :: XNEW(N)
              REAL(KIND=8), INTENT(INOUT) :: XPT(NPT,N)
              REAL(KIND=8), INTENT(INOUT) :: FVAL(NPT)
              REAL(KIND=8), INTENT(INOUT) :: GQ(N)
              REAL(KIND=8), INTENT(INOUT) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(INOUT) :: PQ(NPT)
              REAL(KIND=8), INTENT(INOUT) :: BMAT(NPT+N,N)
              REAL(KIND=8), INTENT(INOUT) :: ZMAT(NPT,NPT-N-1)
              INTEGER(KIND=4), INTENT(IN) :: NDIM
              REAL(KIND=8), INTENT(INOUT) :: D(N)
              REAL(KIND=8), INTENT(INOUT) :: VLAG(NPT+N)
              REAL(KIND=8), INTENT(INOUT) :: W(10*(NPT+N))
              REAL(KIND=8), INTENT(OUT) :: F
              INTEGER(KIND=4), INTENT(OUT) :: INFO
              REAL(KIND=8), INTENT(IN) :: FTARGET
            END SUBROUTINE NEWUOB
          END INTERFACE 
        END MODULE NEWUOB__genmod
