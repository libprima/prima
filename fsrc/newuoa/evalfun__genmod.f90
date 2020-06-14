        !COMPILER-GENERATED INTERFACE MODULE: Mon Jun 15 02:05:56 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EVALFUN__genmod
          INTERFACE 
            SUBROUTINE EVALFUN(N,X,F,XISNAN)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(OUT) :: F
              LOGICAL(KIND=4), INTENT(OUT) :: XISNAN
            END SUBROUTINE EVALFUN
          END INTERFACE 
        END MODULE EVALFUN__genmod
