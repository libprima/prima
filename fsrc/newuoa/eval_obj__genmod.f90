        !COMPILER-GENERATED INTERFACE MODULE: Sun Jun 14 23:59:55 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EVAL_OBJ__genmod
          INTERFACE 
            SUBROUTINE EVAL_OBJ(N,X,F,XISNAN)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(OUT) :: F
              LOGICAL(KIND=4), INTENT(OUT) :: XISNAN
            END SUBROUTINE EVAL_OBJ
          END INTERFACE 
        END MODULE EVAL_OBJ__genmod
