        !COMPILER-GENERATED INTERFACE MODULE: Mon Jun 22 01:39:11 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE HESSMUL__genmod
          INTERFACE 
            SUBROUTINE HESSMUL(N,NPT,XPT,HQ,PQ,D,HD)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: XPT(NPT,N)
              REAL(KIND=8), INTENT(IN) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(IN) :: PQ(NPT)
              REAL(KIND=8), INTENT(IN) :: D(N)
              REAL(KIND=8), INTENT(OUT) :: HD(N)
            END SUBROUTINE HESSMUL
          END INTERFACE 
        END MODULE HESSMUL__genmod
