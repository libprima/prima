        !COMPILER-GENERATED INTERFACE MODULE: Sun Jun 28 01:35:57 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE__genmod
          INTERFACE 
            SUBROUTINE UPDATE(N,NPT,BMAT,ZMAT,IDZ,VLAG,BETA,KNEW)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(INOUT) :: BMAT(N,NPT+N)
              REAL(KIND=8), INTENT(INOUT) :: ZMAT(NPT,NPT-N-1)
              INTEGER(KIND=4), INTENT(INOUT) :: IDZ
              REAL(KIND=8), INTENT(INOUT) :: VLAG(NPT+N)
              REAL(KIND=8), INTENT(IN) :: BETA
              INTEGER(KIND=4), INTENT(IN) :: KNEW
            END SUBROUTINE UPDATE
          END INTERFACE 
        END MODULE UPDATE__genmod
