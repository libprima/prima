! INFNAN is a module defining some procedures concerning INF or NAN.

! Which Fortran standard to follow?
#ifdef __FORTRAN_STANDARD__
#undef __FORTRAN_STANDARD__
#endif
#define __FORTRAN_STANDARD__ 2003 

! Decide whether IEEE_ARITHMETIC is available by checking the version of
! the compiler. 
#ifdef __USE_IEEE_ARITHMETIC__      
#undef __USE_IEEE_ARITHMETIC__
#endif
#define __USE_IEEE_ARITHMETIC__ 0

#if __FORTRAN_STANDARD__ >= 2003 

#ifdef __GFORTRAN__
#if __GNUC__ > 4
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#ifdef __INTEL_COMPILER
#if __INTEL_COMPILER > 1100
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#ifdef __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE > 40
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#ifdef __G95__
#if __G95_MINOR__ > 90
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#endif


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

      pure elemental logical function is_posinf(x)
      implicit none
      real(kind = rp), intent(in) :: x
      is_posinf = (x > hugenum)
      end function is_posinf

      pure elemental logical function is_neginf(x)
      implicit none
      real(kind = rp), intent(in) :: x
      is_neginf = (-x > hugenum)
      end function is_neginf

      pure elemental logical function is_inf(x)
      implicit none
      real(kind = rp), intent(in) :: x
      is_inf = (abs(x) > hugenum)
      end function is_inf

      pure elemental logical function is_nan(x)
      implicit none
      real(kind = rp), intent(in) :: x
#if __USE_IEEE_ARITHMETIC__ == 1
      is_nan = ieee_is_nan(x)
#else
      is_nan = (x /= x)
#endif
      end function is_nan

      pure elemental logical function is_finite(x)
      implicit none
      real(kind = rp), intent(in) :: x
#if __USE_IEEE_ARITHMETIC__ == 1
      is_finite = ieee_is_finite(x)
#else
      is_finite = .not. (is_inf(x) .or. is_nan(x))
#endif
      end function is_finite

      end module infnan
