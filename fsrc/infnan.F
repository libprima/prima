! INFNAN is a module defining some procedures concerning INF or NAN.

#include "ppflags.h"

      module infnan

#if __USE_IEEE_ARITHMETIC__ != 0 
      use ieee_arithmetic, only : is_nan => ieee_is_nan,                &
     & is_finite => ieee_is_finite
#endif

      implicit none
      private
      public :: is_nan, is_finite, is_posinf, is_neginf, is_inf

#if __USE_IEEE_ARITHMETIC__ == 0
      interface is_nan
          module procedure is_nan_32, is_nan_64, is_nan_128
      end interface is_nan
#endif

      contains
      ! The following functions check whether a real number x is
      ! infinite, nan, or finite. Starting from Fortran 2003, the
      ! intrincic ieee_arithmetic module provides ieee_is_nan() and
      ! ieee_is_finite().

#if __USE_IEEE_ARITHMETIC__ == 0 
      pure elemental function is_nan_32(x) result(y)
      use consts, only : REAL32 
      implicit none
      real(REAL32), intent(in) :: x
      logical :: y
      y = (.not. (abs(x) >= 0.0_REAL32))
      ! If we defined is_nan = (x /= x), then gfortran will complain 
      ! about inequality comparison between floating-point numbers.
      end function is_nan_32

      pure elemental function is_nan_64(x) result(y)
      use consts, only : REAL64 
      implicit none
      real(REAL64), intent(in) :: x
      logical :: y
      y = (.not. (abs(x) >= 0.0_REAL64))
      ! If we defined is_nan = (x /= x), then gfortran will complain 
      ! about inequality comparison between floating-point numbers.
      end function is_nan_64

      pure elemental function is_nan_128(x) result(y)
      use consts, only : REAL128 
      implicit none
      real(REAL128), intent(in) :: x
      logical :: y
      y = (.not. (abs(x) >= 0.0_REAL128))
      ! If we defined is_nan = (x /= x), then gfortran will complain 
      ! about inequality comparison between floating-point numbers.
      end function is_nan_128

      pure elemental function is_finite(x) result(y)
      use consts, only : RP, HUGENUM
      implicit none
      real(RP), intent(in) :: x
      logical :: y
      y = (abs(x) <= HUGENUM)
      end function is_finite
#endif

      pure elemental function is_posinf(x) result(y)
      use consts, only : RP, HUGENUM
      implicit none
      real(RP), intent(in) :: x
      logical :: y
      y = (x > HUGENUM)
      end function is_posinf

      pure elemental function is_neginf(x) result(y)
      use consts, only : RP, HUGENUM
      implicit none
      real(RP), intent(in) :: x
      logical :: y
      y = (-x > HUGENUM)
      end function is_neginf

      pure elemental function is_inf(x) result(y)
      use consts, only : RP, HUGENUM
      implicit none
      real(RP), intent(in) :: x
      logical :: y
      y = (abs(x) > HUGENUM)
      end function is_inf

      end module infnan
