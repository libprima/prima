      ! pdfoconst is a module defining some constants to be used by PDFO
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
      
      module pdfoconst 
      ! pdfoconst defines some constants
      
      implicit none
      integer, parameter :: dp = kind(0.0D0)
      ! dp is the kind for double precision
      integer, parameter :: sp = kind(0.0)
      ! sp is the kind for single precision
      real(kind=dp), parameter :: zero = 0.0d0
      real(kind=dp), parameter :: one = 1.0d0
      real(kind=dp), parameter :: two = 2.0d0
      real(kind=dp), parameter :: ten = 10.0d0
      real(kind=dp), parameter :: tenth = 0.1d0
      real(kind=dp), parameter :: half = 0.5d0
      real(kind=dp), parameter :: hugenum = huge(half)
      real(kind=dp), parameter :: hugefun = min(1.0d42, sqrt(hugenum))
      real(kind=dp), parameter :: hugecon = min(1.0d42, sqrt(hugenum))
      
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
