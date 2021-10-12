! test_kind.F
#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      use ieee_arithmetic, only : ieee_is_nan
      implicit none
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

      write(99, *) kind(ieee_is_nan(1.0))  ! This prints 8 in fort.99
      write(99, *) kind(.false.)  ! A comparison, which prints 4 in fort.99
      close(99)

      end subroutine mexFunction
