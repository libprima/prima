! INFNAN is a module defining some procedures concerning INF or NAN.
#include "ppflags.h"

      module infnan

#if __USE_IEEE_ARITHMETIC__ == 1
      use, intrinsic :: ieee_arithmetic, only:ieee_is_nan,ieee_is_finite
#endif

      use consts, only : rp, hugenum
      implicit none

      contains
      ! The following functions check whether a real number x is
      ! infinite, nan, or finite. Starting from Fortran 2003, the
      ! intrincic ieee_arithmetic module provides ieee_is_nan() and
      ! ieee_is_finite().

      pure elemental function is_posinf(x) result(y)
      implicit none
      real(kind = rp), intent(in) :: x
      logical :: y
      y = (x > hugenum)
      end function is_posinf

      pure elemental function is_neginf(x) result(y)
      implicit none
      real(kind = rp), intent(in) :: x
      logical :: y
      y = (-x > hugenum)
      end function is_neginf

      pure elemental function is_inf(x) result(y)
      implicit none
      real(kind = rp), intent(in) :: x
      logical :: y
      y = (abs(x) > hugenum)
      end function is_inf

      pure elemental function is_nan(x) result(y)
      implicit none
      real(kind = rp), intent(in) :: x
      logical :: y
#if __USE_IEEE_ARITHMETIC__ == 1
      y = ieee_is_nan(x)
#else
      y = (.not. (abs(x) >= 0))
      ! Do NOT use is_nan = (x /= x), which will cause gfortran to
      ! complain about equality/inequality comparison between floating
      ! point numbers.
#endif
      end function is_nan

      pure elemental function is_finite(x) result(y)
      implicit none
      real(kind = rp), intent(in) :: x
      logical :: y
#if __USE_IEEE_ARITHMETIC__ == 1
      y = ieee_is_finite(x)
#else
      y = (abs(x) <= hugenum)
#endif
      end function is_finite

      end module infnan
