        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 26 22:34:42 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BIGLAG__genmod
          INTERFACE 
            SUBROUTINE BIGLAG(N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,KNEW,DELTA,D,&
     &ALPHA)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: XOPT(N)
              REAL(KIND=8), INTENT(IN) :: XPT(NPT,N)
              REAL(KIND=8), INTENT(IN) :: BMAT(NPT+N,N)
              REAL(KIND=8), INTENT(IN) :: ZMAT(NPT,NPT-N-1)
              INTEGER(KIND=4), INTENT(IN) :: IDZ
              INTEGER(KIND=4), INTENT(IN) :: KNEW
              REAL(KIND=8), INTENT(IN) :: DELTA
              REAL(KIND=8), INTENT(OUT) :: D(N)
              REAL(KIND=8), INTENT(OUT) :: ALPHA
            END SUBROUTINE BIGLAG
          END INTERFACE 
        END MODULE BIGLAG__genmod
