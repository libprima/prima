        !COMPILER-GENERATED INTERFACE MODULE: Sat Jul 11 01:40:44 2020
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
              REAL(KIND=16), INTENT(IN) :: BMAT(N,NPT+N)
              REAL(KIND=16), INTENT(IN) :: ZMAT(NPT,NPT-N-1)
              REAL(KIND=16), INTENT(IN) :: XPT(N,NPT)
              REAL(KIND=16), INTENT(IN) :: XOPT(N)
              REAL(KIND=16), INTENT(IN) :: D(N)
              REAL(KIND=16), INTENT(OUT) :: VLAG(NPT+N)
              REAL(KIND=16), INTENT(OUT) :: BETA
              REAL(KIND=16), INTENT(OUT) :: WCHECK(NPT)
              REAL(KIND=16) :: DSQ
              REAL(KIND=16) :: XOPTSQ
            END SUBROUTINE VLAGBETA
          END INTERFACE 
        END MODULE VLAGBETA__genmod
