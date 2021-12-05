        !COMPILER-GENERATED INTERFACE MODULE: Sat Dec  4 12:27:35 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INIT__genmod
          INTERFACE 
            SUBROUTINE INIT(PROB,PROBNAME,N)
              USE TEST
              TYPE (PROBLEM), INTENT(INOUT) :: PROB
              CHARACTER(*), INTENT(IN) :: PROBNAME
              INTEGER(KIND=4), INTENT(IN) :: N
            END SUBROUTINE INIT
          END INTERFACE 
        END MODULE INIT__genmod
