        !COMPILER-GENERATED INTERFACE MODULE: Sun Jun 28 00:20:14 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATEQ__genmod
          INTERFACE 
            SUBROUTINE UPDATEQ(N,NPT,IDZ,KNEW,FQDIFF,XPTKNEW,BMATKNEW,  &
     &ZMAT,GQ,HQ,PQ)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: IDZ
              INTEGER(KIND=4), INTENT(IN) :: KNEW
              REAL(KIND=8), INTENT(IN) :: FQDIFF
              REAL(KIND=8), INTENT(IN) :: XPTKNEW(N)
              REAL(KIND=8), INTENT(IN) :: BMATKNEW(N)
              REAL(KIND=8), INTENT(IN) :: ZMAT(NPT,NPT-N-1)
              REAL(KIND=8), INTENT(INOUT) :: GQ(N)
              REAL(KIND=8), INTENT(INOUT) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(INOUT) :: PQ(NPT)
            END SUBROUTINE UPDATEQ
          END INTERFACE 
        END MODULE UPDATEQ__genmod
