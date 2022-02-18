#include "fintrf.h"
C======================================================================
#if 0
C     
C     timestwo.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     timestwo.f
C
C     Computational function that takes a scalar and doubles it.
      
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2018 The MathWorks, Inc.
C     
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:

#if MX_HAS_INTERLEAVED_COMPLEX
      mwPointer mxGetDoubles
#else
      mwPointer mxGetPr
#endif

      mwPointer mxCreateDoubleMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN

C     Pointers to input/output mxArrays:
      mwPointer x_ptr, y_ptr

C     Array information:
      mwPointer mrows, ncols
      mwSize size

C     Arguments for computational routine:
      real*8  x_input, y_output

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:nInput',
     +                           'One input required.')
      elseif(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:nOutput',
     +                           'Too many output arguments.')
      endif

C     Validate inputs
C     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:timestwo:NonNumeric',
     +                           'Input must be a number.')
      endif

C     Get the size of the input array.
      mrows = mxGetM(prhs(1))
      ncols = mxGetN(prhs(1))
      size = mrows*ncols

C     Create Fortran array from the input argument.

#if MX_HAS_INTERLEAVED_COMPLEX
      x_ptr = mxGetDoubles(prhs(1))
#else
      x_ptr = mxGetPr(prhs(1))
#endif
      call mxCopyPtrToReal8(x_ptr,x_input,size)

C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(mrows,ncols,0)
#if MX_HAS_INTERLEAVED_COMPLEX
      y_ptr = mxGetDoubles(plhs(1))
#else
      y_ptr = mxGetPr(plhs(1))
#endif

C     Call the computational subroutine.
      call timestwo(y_output, x_input)

C     Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyReal8ToPtr(y_output,y_ptr,size)     

      return
      end

C-----------------------------------------------------------------------
C     Computational routine

      subroutine timestwo(y_output, x_input)
      real*8 x_input, y_output

      y_output = 2.0 * x_input
      return
      end
