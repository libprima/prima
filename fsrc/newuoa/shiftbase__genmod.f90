        !COMPILER-GENERATED INTERFACE MODULE: Sun Jun 28 00:20:13 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SHIFTBASE__genmod
          INTERFACE 
            SUBROUTINE SHIFTBASE(N,NPT,IDZ,XOPT,PQ,BMAT,ZMAT,GQ,HQ,XPT)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: IDZ
              REAL(KIND=8), INTENT(IN) :: XOPT(N)
              REAL(KIND=8), INTENT(IN) :: PQ(NPT)
              REAL(KIND=8), INTENT(INOUT) :: BMAT(N,NPT+N)
              REAL(KIND=8), INTENT(IN) :: ZMAT(NPT,NPT-N-1)
              REAL(KIND=8), INTENT(INOUT) :: GQ(N)
              REAL(KIND=8), INTENT(INOUT) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(INOUT) :: XPT(N,NPT)
            END SUBROUTINE SHIFTBASE
          END INTERFACE 
        END MODULE SHIFTBASE__genmod
