        !COMPILER-GENERATED INTERFACE MODULE: Sat Jul 11 01:40:43 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NEWUOA__genmod
          INTERFACE 
            SUBROUTINE NEWUOA(NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,F,INFO, &
     &FTARGET)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              REAL(KIND=16), INTENT(INOUT) :: X(:)
              REAL(KIND=16), INTENT(IN) :: RHOBEG
              REAL(KIND=16), INTENT(IN) :: RHOEND
              INTEGER(KIND=4), INTENT(IN) :: IPRINT
              INTEGER(KIND=4), INTENT(IN) :: MAXFUN
              REAL(KIND=16), INTENT(OUT) :: F
              INTEGER(KIND=4), INTENT(OUT) :: INFO
              REAL(KIND=16), INTENT(IN) :: FTARGET
            END SUBROUTINE NEWUOA
          END INTERFACE 
        END MODULE NEWUOA__genmod
