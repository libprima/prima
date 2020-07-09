#include "fintrf.h"

      module mex_mod

      use consts, only : INT32, DP, RP
      implicit none
      private
      public :: notComplex, mwOne, intOne, intTwo, convTol
      public :: mexErrMsgIdAndTxt
      public :: mxCopyReal8ToPtr
      public :: mxCopyPtrToReal8
      public :: mxDestroyArray
      public :: mexCallMATLAB
      public :: mxCreateDoubleMatrix
      public :: mxCreateDoubleScalar
      public :: mxGetDoubles
      public :: mxGetPr
      public :: mxGetM
      public :: mxGetN
      public :: mxIsDouble
      public :: mxIsClass

      
      integer(INT32), parameter :: notComplex = 0
C notComplex is used in mxCreateDoubleMatrix
      integer(INT32), parameter :: intOne=1, intTwo=2 
C intOne, intTwo, funFeval are  used when calling mexCallMATLAB 
      mwSize, parameter :: mwOne = 1 ! Integer 1 with type mwSize
C mwOne is used in mxCreateDoubleMatrix and mxCopyPtrToReal8
      real(DP), parameter :: convTol = 1.0e1_DP*max(epsilon(0.0_DP),
     &    real(epsilon(0.0_RP), DP))
C convTol is the tolerance of difference due to real type conversion


      interface 

          subroutine mexErrMsgIdAndTxt(errorid, errormsg)
              implicit none
              character*(*), intent(in) ::  errorid, errormsg
          end subroutine mexErrMsgIdAndTxt
    
          subroutine  mxCopyPtrToReal8(px, y, n)
              use consts, only : DP
              implicit none
              mwPointer, intent(in) :: px
              mwSize, intent(in) :: n
              real(DP), intent(out) :: y(n)
          end subroutine mxCopyPtrToReal8
    
          subroutine mxCopyReal8ToPtr(y, px, n)
              use consts, only : DP
              implicit none
              mwPointer, intent(in) :: px
              mwSize, intent(in) :: n
              real(DP), intent(out) :: y(n)
          end subroutine mxCopyReal8ToPtr
    
          subroutine mxDestroyArray(pm)
              implicit none
              mwPointer, intent(in) :: pm
          end subroutine mxDestroyArray
    
    
          function mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName)
              use consts, only : INT32
              implicit none
              integer(INT32) :: mexCallMATLAB
              integer(INT32), intent(in) :: nlhs, nrhs
              mwPointer, intent(in) :: plhs(*), prhs(*)
              character*(*), intent(in) :: functionName
          end function mexCallMATLAB
    
          function mxCreateDoubleMatrix(m, n, ComplexFlag)
              use consts, only : INT32
              implicit none
              mwPointer :: mxCreateDoubleMatrix
              mwSize, intent(in) :: m, n
              integer(INT32), intent(in) :: ComplexFlag
          end function mxCreateDoubleMatrix
    
          function mxCreateDoubleScalar(x)
              use consts, only : DP
              implicit none
              mwPointer :: mxCreateDoubleScalar
              real(DP), intent(in) :: x
          end function mxCreateDoubleScalar
    
          function mxGetDoubles(pa)
              implicit none
              mwPointer :: mxGetDoubles                  
              mwPointer, intent(in) :: pa
          end function mxGetDoubles
    
          function mxGetPr(pa)
              implicit none
              mwPointer :: mxGetPr                  
              mwPointer, intent(in) :: pa
          end function mxGetPr
    
          function mxGetM(pm)
              implicit none
              mwPointer :: mxGetM
              mwPointer, intent(in) :: pm
          end function mxGetM
    
          function mxGetN(pm)
              implicit none
              mwPointer :: mxGetN
              mwPointer, intent(in) :: pm
          end function mxGetN
    
          function mxIsDouble(pm)
              use consts, only : INT32
              implicit none
              integer(INT32) :: mxIsDouble
              mwPointer, intent(in) :: pm
          end function mxIsDouble
    
          function mxIsClass(pm, classname)
              use consts, only : INT32
              implicit none
              integer(INT32) :: mxIsClass
              mwPointer, intent(in) :: pm
              character*(*), intent(in) :: classname
          end function mxIsClass

      end interface

      end module mex_mod
