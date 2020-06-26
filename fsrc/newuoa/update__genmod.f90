        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun 27 00:16:41 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE__genmod
          INTERFACE 
            SUBROUTINE UPDATE(N,NPT,BMAT,ZMAT,IDZ,VLAG,BETA,KNEW)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(INOUT) :: BMAT(NPT+N,N)
              REAL(KIND=8), INTENT(INOUT) :: ZMAT(NPT,NPT-N-1)
              INTEGER(KIND=4), INTENT(INOUT) :: IDZ
              REAL(KIND=8), INTENT(INOUT) :: VLAG(NPT+N)
              REAL(KIND=8), INTENT(IN) :: BETA
              INTEGER(KIND=4), INTENT(IN) :: KNEW
            END SUBROUTINE UPDATE
          END INTERFACE 
        END MODULE UPDATE__genmod
