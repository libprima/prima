        !COMPILER-GENERATED INTERFACE MODULE: Sat Jul 11 11:20:37 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MEXFUNCTION__genmod
          INTERFACE 
            SUBROUTINE MEXFUNCTION(NLHS,PLHS,NRHS,PRHS)
              INTEGER(KIND=4), INTENT(IN) :: NRHS
              INTEGER(KIND=4), INTENT(IN) :: NLHS
              INTEGER(KIND=8), INTENT(OUT) :: PLHS(NLHS)
              INTEGER(KIND=8), INTENT(IN) :: PRHS(NRHS)
            END SUBROUTINE MEXFUNCTION
          END INTERFACE 
        END MODULE MEXFUNCTION__genmod
