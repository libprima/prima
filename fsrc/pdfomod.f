      ! pdfomod is a module defining some constants and
      ! functions/subroutines to be used by PDFO
      !
      !*****************************************************************
      !   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk) 
      !               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
      !               Department of Applied Mathematics,
      !               The Hong Kong Polytechnic University
      !
      !   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
      !****************************************************************

      ! Remarks:
      !
      ! 1. REAL*4, REAL*8, INTEGER*4, INTEGER*8 are not Fortran standard
      !    expressions. Although they are supported by many compilers as
      !    extensions, it is better to avoid them.
      !
      ! 2. Never use KIND with a literal value (e.g., REAL(KIND=8)),
      !    because Fortran standards never define what KIND=8 means.
      !    There is NO guarantee that REAL(KIND=8) will be of double
      !    precision.

      module pdfomod 
      ! pdfomod defines some constants and functions/subroutines

      implicit none
      integer, parameter :: dp = kind(0.0d0), sp = kind(0.0)
      ! dp is the kind for double precision
      ! sp is the kind for single precision
      integer, parameter :: rp = dp  
      ! rp is the kind for the default real precision
      real(kind=rp), parameter :: zero = 0.0_rp
      real(kind=rp), parameter :: one = 1.0_rp
      real(kind=rp), parameter :: two = 2.0_rp
      real(kind=rp), parameter :: half = 0.5_rp
      real(kind=rp), parameter :: ten = 10.0_rp
      real(kind=rp), parameter :: tenth = 0.1_rp
      real(kind=rp), parameter :: pi = real(acos(-1.0_dp), rp)
      real(kind=rp), parameter :: hugenum = huge(zero)
      real(kind=rp), parameter :: almost_infinity = half*hugenum 
      real(kind=rp), parameter :: hugefun = min(1.0e42_rp,sqrt(hugenum))
      real(kind=rp), parameter :: hugecon = min(1.0e42_rp,sqrt(hugenum))

      integer, parameter :: int4 = selected_int_kind(8)
      ! int4 is the kind for integer*4
      ! SELECTED_INT_KIND(p) returns an INTEGER that equals the
      ! processor dependent kind type parameter of the integer type 
      ! accommodating all values n with -10^p < n < 10^p. Therefore,
      ! SELECTED_INT_KIND(p) should be the kind of integer*4 for p = 7,
      ! 8, 9 (also for p = 5 and 6 unless integer*3 is supportred).

      contains

      ! The following functions check whether a real number x is
      ! infinite, nan, or finite. Starting from Fortran 2003, the
      ! intrincic ieee_arithmetic module provides ieee_is_nan() and
      ! ieee_is_finite(); gfortran provides isnan() as an extension.
      pure elemental logical function is_inf(x)
      implicit none
      real(kind = rp), intent(in) :: x
      is_inf = (x > hugenum)
      end function is_inf

      pure elemental logical function is_nan(x)
      implicit none
      real(kind = rp), intent(in) :: x
      is_nan = (x /= x)
      end function is_nan

      pure elemental logical function is_finite(x)
      implicit none
      real(kind = rp), intent(in) :: x
      is_finite = .not. (is_inf(x) .or. is_nan(x))
      end function is_finite

      end module pdfomod
