        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 25 14:58:28 2021
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRSTLP__genmod
          INTERFACE 
            SUBROUTINE TRSTLP(N,M,A,B,RHO,DX,IFULL,IACT,Z,ZDOTA,VMULTC, &
     &SDIRN,DXNEW,VMULTD)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: M
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(IN) :: RHO
              REAL(KIND=8), INTENT(INOUT) :: DX(:)
              INTEGER(KIND=4), INTENT(OUT) :: IFULL
              INTEGER(KIND=4), INTENT(INOUT) :: IACT(:)
              REAL(KIND=8), INTENT(INOUT) :: Z(:,:)
              REAL(KIND=8), INTENT(INOUT) :: ZDOTA(:)
              REAL(KIND=8), INTENT(INOUT) :: VMULTC(:)
              REAL(KIND=8), INTENT(INOUT) :: SDIRN(:)
              REAL(KIND=8), INTENT(INOUT) :: DXNEW(:)
              REAL(KIND=8), INTENT(INOUT) :: VMULTD(:)
            END SUBROUTINE TRSTLP
          END INTERFACE 
        END MODULE TRSTLP__genmod
