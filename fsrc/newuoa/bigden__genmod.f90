        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 24 03:07:14 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BIGDEN__genmod
          INTERFACE 
            SUBROUTINE BIGDEN(N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,KOPT,KNEW,D, &
     &WCHECK,VLAG,BETA)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: XOPT(N)
              REAL(KIND=8), INTENT(IN) :: XPT(NPT,N)
              REAL(KIND=8), INTENT(IN) :: BMAT(NPT+N,N)
              REAL(KIND=8), INTENT(IN) :: ZMAT(NPT,NPT-N-1)
              INTEGER(KIND=4), INTENT(IN) :: IDZ
              INTEGER(KIND=4), INTENT(IN) :: KOPT
              INTEGER(KIND=4), INTENT(IN) :: KNEW
              REAL(KIND=8), INTENT(INOUT) :: D(N)
              REAL(KIND=8), INTENT(OUT) :: WCHECK(NPT+N)
              REAL(KIND=8), INTENT(OUT) :: VLAG(NPT+N)
              REAL(KIND=8), INTENT(OUT) :: BETA
            END SUBROUTINE BIGDEN
          END INTERFACE 
        END MODULE BIGDEN__genmod
