        !COMPILER-GENERATED INTERFACE MODULE: Sat Jul 11 01:40:44 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALQUAD__genmod
          INTERFACE 
            SUBROUTINE CALQUAD(VQUAD,D,X,XPT,GQ,HQ,PQ,N,NPT,WCHECK)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=16), INTENT(OUT) :: VQUAD
              REAL(KIND=16), INTENT(IN) :: D(N)
              REAL(KIND=16), INTENT(IN) :: X(N)
              REAL(KIND=16), INTENT(IN) :: XPT(N,NPT)
              REAL(KIND=16), INTENT(IN) :: GQ(N)
              REAL(KIND=16), INTENT(IN) :: HQ(N,N)
              REAL(KIND=16), INTENT(IN) :: PQ(NPT)
              REAL(KIND=16) :: WCHECK(NPT)
            END SUBROUTINE CALQUAD
          END INTERFACE 
        END MODULE CALQUAD__genmod
