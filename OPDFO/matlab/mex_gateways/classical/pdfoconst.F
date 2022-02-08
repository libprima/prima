C pdfoconst is a module that defines some constants to be used by PDFO
C
C **********************************************************************
C   Authors:    Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk) 
C               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
C               Department of Applied Mathematics,
C               The Hong Kong Polytechnic University
C
C   Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
C **********************************************************************

C Remarks:
C
C 1. REAL*4, REAL*8, INTEGER*4, INTEGER*8 are not Fortran standard 
C expressions. Although they are supported by many compilers as extensions, 
C it is better to avoid them.
C
C 2. Never use KIND with a literal value (e.g., REAL(KIND=8)), because 
C Fortran standards never define what KIND=8 means. There is NO guarantee 
C that REAL(KIND=8) will be of double precision.

      module pdfoconst 
      ! pdfoconst defines some constants

      implicit none
      integer, parameter :: dp = kind(0.0D0)
      ! dp is the kind for double precision
      integer, parameter :: sp = kind(0.0)
      ! sp is the kind for single precision
      real(kind=dp), parameter :: ZERO = 0.0D0
      real(kind=dp), parameter :: ONE = 1.0D0
      real(kind=dp), parameter :: TWO = 2.0D0
      real(kind=dp), parameter :: TEN = 10.0D0
      real(kind=dp), parameter :: TENTH = 0.1D0
      real(kind=dp), parameter :: HALF = 0.5D0
      real(kind=dp), parameter :: HUGENUM = huge(HALF)
      integer, parameter :: maxe = maxexponent(zero)
      real(kind=dp), parameter :: hugefun = 2.0D0**min(100, maxe/2)
      real(kind=dp), parameter :: hugecon = hugefun

      integer, parameter :: int4 = selected_int_kind(8)
      ! int4 is the kind for integer*4
      ! SELECTED_INT_KIND(p) returns an INTEGER that equals the
      ! processor dependent kind type parameter of the integer type 
      ! accommodating all values n with -10^p < n < 10^p. Therefore,
      ! SELECTED_INT_KIND(p) should be the kind of integer*4 for p = 7,
      ! 8, 9 (also for p = 5 and 6 unless integer*3 is supportred).
      
      end module pdfoconst

      subroutine use_pdfoconst ()
      ! This is a function that does nothing. It is to entertain F2PY,
      ! which interprets a file only if it contains at least one
      ! function or subroutine (as of NumPy v1.17).

      use pdfoconst

      return
      end subroutine use_pdfoconst
