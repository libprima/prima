#include "fintrf.h"

      module mex_mod

      use consts_mod, only : INT32, DP, RP
      implicit none
      private

      public :: notComplex, mwOne, intOne, intTwo, convTol

      ! MEX subroutines
      public :: mexErrMsgIdAndTxt
      public :: mxCopyPtrToReal8
      public :: mxCopyReal8ToPtr
      public :: mxDestroyArray

      ! MEX functions 
      public :: mexCallMATLAB
      public :: mxCreateDoubleMatrix
      public :: mxCreateDoubleScalar
      public :: mxGetDoubles
      public :: mxGetM
      public :: mxGetN
      public :: mxGetPr
      public :: mxGetString
      public :: mxIsChar
      public :: mxIsClass
      public :: mxIsDouble

      
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
              use consts_mod, only : DP
              implicit none
              mwPointer, intent(in) :: px
              mwSize, intent(in) :: n
              real(DP), intent(out) :: y(n)
          end subroutine mxCopyPtrToReal8
    
          subroutine mxCopyReal8ToPtr(y, px, n)
              use consts_mod, only : DP
              implicit none
              mwPointer, intent(in) :: px
              mwSize, intent(in) :: n
              real(DP), intent(in) :: y(n)
          end subroutine mxCopyReal8ToPtr
    
          subroutine mxDestroyArray(pm)
              implicit none
              mwPointer, intent(in) :: pm
          end subroutine mxDestroyArray
    
    
          function mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName)
              use consts_mod, only : INT32
              implicit none
              integer(INT32) :: mexCallMATLAB
              integer(INT32), intent(in) :: nlhs, nrhs
              ! N.B.: Segmentation Fault will occur if we write 
              !mwPointer, intent(in) :: plhs(:), prhs(:) 
              mwPointer, intent(in) :: plhs(*), prhs(*) 
              character*(*), intent(in) :: functionName
          end function mexCallMATLAB
    
          function mxCreateDoubleMatrix(m, n, ComplexFlag)
              use consts_mod, only : INT32
              implicit none
              mwPointer :: mxCreateDoubleMatrix
              mwSize, intent(in) :: m, n
              integer(INT32), intent(in) :: ComplexFlag
          end function mxCreateDoubleMatrix
    
          function mxCreateDoubleScalar(x)
              use consts_mod, only : DP
              implicit none
              mwPointer :: mxCreateDoubleScalar
              real(DP), intent(in) :: x
          end function mxCreateDoubleScalar
    
          function mxGetDoubles(pa)
              implicit none
              mwPointer :: mxGetDoubles                  
              mwPointer, intent(in) :: pa
          end function mxGetDoubles
    
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

          function mxGetPr(pa)
              implicit none
              mwPointer :: mxGetPr                  
              mwPointer, intent(in) :: pa
          end function mxGetPr

          function mxGetString(pm, str, strlen)
              use consts_mod, only : INT32
              implicit none
              integer(INT32) :: mxGetString
              mwPointer, intent(in) :: pm
              character*(*), intent(in) :: str
              mwSize, intent(in) :: strlen
          end function mxGetString

          function mxIsClass(pm, classname)
              use consts_mod, only : INT32
              implicit none
              integer(INT32) :: mxIsClass
              mwPointer, intent(in) :: pm
              character*(*), intent(in) :: classname
          end function mxIsClass

          function mxIsChar(pm)
              use consts_mod, only : INT32
              implicit none
              integer(INT32) :: mxIsChar
              mwPointer, intent(in) :: pm
          end function mxIsChar

          function mxIsDouble(pm)
              use consts_mod, only : INT32
              implicit none
              integer(INT32) :: mxIsDouble
              mwPointer, intent(in) :: pm
          end function mxIsDouble
    
      end interface

      end module mex_mod

      ! Remark: What is the distinction between mx and mex prefixes?
      ! According to 
      ! matlab.izmiran.ru/help/techdoc/matlab_external/ch03cre5.html
      ! "Routines in the API that are prefixed with mx allow you to 
      ! create, access, manipulate, and destroy mxArrays. Routines 
      ! prefixed with mex perform operations back in the MATLAB 
      ! environment."
