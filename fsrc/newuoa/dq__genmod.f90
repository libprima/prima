        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 16 22:22:46 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DQ__genmod
          INTERFACE 
            SUBROUTINE DQ(QDIFF,D,X,XPT,GQ,HQ,PQ,N,NPT)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(OUT) :: QDIFF
              REAL(KIND=8), INTENT(IN) :: D(N)
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(IN) :: XPT(NPT,N)
              REAL(KIND=8), INTENT(IN) :: GQ(N)
              REAL(KIND=8), INTENT(IN) :: HQ(N*(N+1)/2)
              REAL(KIND=8), INTENT(IN) :: PQ(NPT)
            END SUBROUTINE DQ
          END INTERFACE 
        END MODULE DQ__genmod
