        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun 27 00:16:39 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALQUAD__genmod
          INTERFACE 
            SUBROUTINE CALQUAD(VQUAD,D,X,XPT,GQ,HQ,PQ,N,NPT,WCHECK)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(OUT) :: VQUAD
              REAL(KIND=8), INTENT(IN) :: D(N)
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(IN) :: XPT(NPT,N)
              REAL(KIND=8), INTENT(IN) :: GQ(N)
              REAL(KIND=8), INTENT(IN) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(IN) :: PQ(NPT)
              REAL(KIND=8) :: WCHECK(NPT)
            END SUBROUTINE CALQUAD
          END INTERFACE 
        END MODULE CALQUAD__genmod
