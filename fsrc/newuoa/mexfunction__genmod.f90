        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 11 00:31:09 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MEXFUNCTION__genmod
          INTERFACE 
            SUBROUTINE MEXFUNCTION(NARGOUT,POUTPUT,NARGIN,PINPUT)
              INTEGER(KIND=4), INTENT(IN) :: NARGIN
              INTEGER(KIND=4), INTENT(IN) :: NARGOUT
              INTEGER(KIND=8), INTENT(OUT) :: POUTPUT(NARGOUT)
              INTEGER(KIND=8), INTENT(IN) :: PINPUT(NARGIN)
            END SUBROUTINE MEXFUNCTION
          END INTERFACE 
        END MODULE MEXFUNCTION__genmod
