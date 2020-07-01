        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul  1 15:19:13 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BIGLAG__genmod
          INTERFACE 
            SUBROUTINE BIGLAG(X,XPT,BMAT,ZMAT,IDZ,KNEW,DELTA,D,ALPHA)
              REAL(KIND=8), INTENT(IN) :: X(SIZE(XPT,1))
              REAL(KIND=8), INTENT(IN) :: XPT(:,:)
              REAL(KIND=8), INTENT(IN) :: BMAT(SIZE(X),SIZE(X)+SIZE(XPT,&
     &2))
              REAL(KIND=8), INTENT(IN) :: ZMAT(SIZE(XPT,2),SIZE(XPT,2)- &
     &SIZE(X)-1)
              INTEGER(KIND=4), INTENT(IN) :: IDZ
              INTEGER(KIND=4), INTENT(IN) :: KNEW
              REAL(KIND=8), INTENT(IN) :: DELTA
              REAL(KIND=8), INTENT(OUT) :: D(SIZE(X))
              REAL(KIND=8), INTENT(OUT) :: ALPHA
            END SUBROUTINE BIGLAG
          END INTERFACE 
        END MODULE BIGLAG__genmod
