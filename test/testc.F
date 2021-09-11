! testc.F
#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      use ieee_arithmetic, only : ieee_is_nan
      use iso_fortran_env, only : real64, real128
      implicit none
      real(real64) :: a
      real(real128) :: b
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

      a = real(0.59843577329095299D0, real64)
      write(99,*) a
      b = 0.59843577329095299_real128

      write(99, *) 'cos(a) = ',
     & cos(a)
      write(99, *) 'cos(b) = ',
     & cos(b)

      write(99, *) 'cos(0.59843577329095299_real64) = ',
     & cos(0.59843577329095299_real64)
      write(99, *) 'cos(0.59843577329095299_real128) = ',
     & cos(0.59843577329095299_real128)
      close(99)

      end subroutine mexFunction
