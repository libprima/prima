        !COMPILER-GENERATED INTERFACE MODULE: Sat Aug  1 09:20:27 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MEXFUNCTION__genmod
          INTERFACE 
            SUBROUTINE MEXFUNCTION(NLHS,PLHS,NRHS,PRHS)
              INTEGER(KIND=4), INTENT(IN) :: NRHS
              INTEGER(KIND=4), INTENT(IN) :: NLHS
              INTEGER(KIND=8), INTENT(INOUT) :: PLHS(NLHS)
              INTEGER(KIND=8), INTENT(IN) :: PRHS(NRHS)
            END SUBROUTINE MEXFUNCTION
          END INTERFACE 
        END MODULE MEXFUNCTION__genmod
