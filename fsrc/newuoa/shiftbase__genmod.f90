        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun 20 19:12:39 2020
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
              REAL(KIND=8), INTENT(INOUT) :: BMAT(NPT+N,N)
              REAL(KIND=8), INTENT(INOUT) :: ZMAT(NPT,NPT-N-1)
              REAL(KIND=8), INTENT(INOUT) :: GQ(N)
              REAL(KIND=8), INTENT(INOUT) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(INOUT) :: XPT(NPT,N)
            END SUBROUTINE SHIFTBASE
          END INTERFACE 
        END MODULE SHIFTBASE__genmod
