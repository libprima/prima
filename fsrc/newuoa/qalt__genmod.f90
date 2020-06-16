        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 17 01:59:21 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE QALT__genmod
          INTERFACE 
            SUBROUTINE QALT(GQ,HQ,PQ,FVAL,SMAT,ZMAT,N,NPT,KOPT,IDZ)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(OUT) :: GQ(N)
              REAL(KIND=8), INTENT(OUT) :: HQ((N*(N+1))/2)
              REAL(KIND=8), INTENT(OUT) :: PQ(NPT)
              REAL(KIND=8), INTENT(IN) :: FVAL(NPT)
              REAL(KIND=8), INTENT(IN) :: SMAT(NPT,N)
              REAL(KIND=8), INTENT(IN) :: ZMAT(NPT,NPT-N-1)
              INTEGER(KIND=4), INTENT(IN) :: KOPT
              INTEGER(KIND=4), INTENT(IN) :: IDZ
            END SUBROUTINE QALT
          END INTERFACE 
        END MODULE QALT__genmod
