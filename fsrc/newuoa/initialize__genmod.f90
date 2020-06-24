        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 25 01:03:17 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALIZE__genmod
          INTERFACE 
            SUBROUTINE INITIALIZE(N,NPT,RHOBEG,X,XBASE,XPT,F,FVAL,XOPT, &
     &FOPT,KOPT,BMAT,ZMAT,GQ,HQ,PQ,NF,INFO,FTARGET)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: RHOBEG
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(OUT) :: XBASE(N)
              REAL(KIND=8), INTENT(OUT) :: XPT(NPT,N)
              REAL(KIND=8), INTENT(OUT) :: F
              REAL(KIND=8), INTENT(OUT) :: FVAL(NPT)
              REAL(KIND=8), INTENT(OUT) :: XOPT(N)
              REAL(KIND=8), INTENT(OUT) :: FOPT
              INTEGER(KIND=4), INTENT(OUT) :: KOPT
              REAL(KIND=8), INTENT(OUT) :: BMAT(NPT+N,N)
              REAL(KIND=8), INTENT(OUT) :: ZMAT(NPT,NPT-N-1)
              REAL(KIND=8), INTENT(OUT) :: GQ(N)
              REAL(KIND=8), INTENT(OUT) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(OUT) :: PQ(NPT)
              INTEGER(KIND=4), INTENT(OUT) :: NF
              INTEGER(KIND=4), INTENT(OUT) :: INFO
              REAL(KIND=8), INTENT(IN) :: FTARGET
            END SUBROUTINE INITIALIZE
          END INTERFACE 
        END MODULE INITIALIZE__genmod
