        !COMPILER-GENERATED INTERFACE MODULE: Sun Jun 28 00:20:14 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRSAPP__genmod
          INTERFACE 
            SUBROUTINE TRSAPP(N,NPT,TOL,X,XPT,GQ,HQ,PQ,DELTA,S,CRVMIN,  &
     &QRED,INFO)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(IN) :: XPT(N,NPT)
              REAL(KIND=8), INTENT(IN) :: GQ(N)
              REAL(KIND=8), INTENT(IN) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(IN) :: PQ(NPT)
              REAL(KIND=8), INTENT(IN) :: DELTA
              REAL(KIND=8), INTENT(OUT) :: S(N)
              REAL(KIND=8), INTENT(OUT) :: CRVMIN
              REAL(KIND=8), INTENT(OUT) :: QRED
              INTEGER(KIND=4), INTENT(OUT) :: INFO
            END SUBROUTINE TRSAPP
          END INTERFACE 
        END MODULE TRSAPP__genmod
