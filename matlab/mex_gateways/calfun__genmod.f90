        !COMPILER-GENERATED INTERFACE MODULE: Sat Jul 11 11:20:37 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALFUN__genmod
          INTERFACE 
            SUBROUTINE CALFUN(N,X,FUNVAL)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(OUT) :: FUNVAL
            END SUBROUTINE CALFUN
          END INTERFACE 
        END MODULE CALFUN__genmod
