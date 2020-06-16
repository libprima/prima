        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 17 01:59:21 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VLAGBETA__genmod
          INTERFACE 
            SUBROUTINE VLAGBETA(N,NPT,IDZ,KOPT,BMAT,ZMAT,XPT,XOPT,D,VLAG&
     &,BETA,WCHECK,DSQ,XOPTSQ)
              INTEGER(KIND=4), INTENT(IN) :: NPT
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: IDZ
              INTEGER(KIND=4), INTENT(IN) :: KOPT
              REAL(KIND=8), INTENT(IN) :: BMAT(NPT+N,N)
              REAL(KIND=8), INTENT(IN) :: ZMAT(NPT,NPT-N-1)
              REAL(KIND=8), INTENT(IN) :: XPT(NPT,N)
              REAL(KIND=8), INTENT(IN) :: XOPT(N)
              REAL(KIND=8), INTENT(IN) :: D(N)
              REAL(KIND=8), INTENT(OUT) :: VLAG(NPT+N)
              REAL(KIND=8), INTENT(OUT) :: BETA
              REAL(KIND=8), INTENT(OUT) :: WCHECK(NPT)
              REAL(KIND=8) :: DSQ
              REAL(KIND=8) :: XOPTSQ
            END SUBROUTINE VLAGBETA
          END INTERFACE 
        END MODULE VLAGBETA__genmod
