!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of trstlp__genmod.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 01-Jul-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!COMPILER-GENERATED INTERFACE MODULE: Fri Jun 25 14:58:28 2021
! This source file is for reference only and may not completely
! represent the generated interface used by the compiler.
              MODULE TRSTLP__genmod
                INTERFACE
                  SUBROUTINE TRSTLP(N,M,A,B,RHO,DX,IFULL,IACT,Z,ZDOTA,VM&
     &ULTC, SDIRN,DXNEW,VMULTD)
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