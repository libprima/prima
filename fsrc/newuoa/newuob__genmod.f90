        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 25 01:48:17 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NEWUOB__genmod
          INTERFACE 
            SUBROUTINE NEWUOB(N,NPT,RHOBEG,RHOEND,IPRINT,MAXFUN,FTARGET,&
     &X,F,NF,INFO)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: RHOBEG
              REAL(KIND=8), INTENT(IN) :: RHOEND
              INTEGER(KIND=4), INTENT(IN) :: IPRINT
              INTEGER(KIND=4), INTENT(IN) :: MAXFUN
              REAL(KIND=8), INTENT(IN) :: FTARGET
              REAL(KIND=8), INTENT(INOUT) :: X(N)
              REAL(KIND=8), INTENT(OUT) :: F
              INTEGER(KIND=4), INTENT(OUT) :: NF
              INTEGER(KIND=4), INTENT(OUT) :: INFO
            END SUBROUTINE NEWUOB
          END INTERFACE 
        END MODULE NEWUOB__genmod
