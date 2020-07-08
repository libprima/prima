! CONSTS is a module defining some constants 

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

#include "ppflags.h"

      module consts 

      implicit none
      private
      public :: DEBUG_MODE
      public :: DP, SP, RP, INT32
      public :: ZERO, ONE, TWO, HALF, QUART, TEN, TENTH, PI
      public :: EPS, HUGENUM, ALMOST_INFINITY, HUGEFUN, HUGECON 

#if __DEBUG_MODE__ == 1
      logical, parameter :: DEBUG_MODE = .true.
#else
      logical, parameter :: DEBUG_MODE = .false.
#endif

      integer, parameter :: DP = kind(0.0d0), SP = kind(0.0)
      ! DP is the kind for double precision
      ! SP is the kind for single precision
      integer, parameter :: RP = DP
      ! RP is the kind for the real precision that will be used

      real(RP), parameter :: ZERO = 0.0_RP
      real(RP), parameter :: ONE = 1.0_RP
      real(RP), parameter :: TWO = 2.0_RP
      real(RP), parameter :: HALF = 0.5_RP
      real(RP), parameter :: QUART = 0.25_RP
      real(RP), parameter :: TEN = 10.0_RP
      real(RP), parameter :: TENTH = 0.1_RP
      real(RP), parameter :: PI = 3.141592653589793238462643383_RP

      real(RP), parameter :: EPS = epsilon(ZERO)
      real(RP), parameter :: HUGENUM = huge(ZERO)
      real(RP), parameter :: ALMOST_INFINITY = HALF*HUGENUM
      
      integer, parameter :: MAXE = maxexponent(ZERO)
      real(RP), parameter :: HUGEFUN = TWO**min(100, MAXE/2)
      real(RP), parameter :: HUGECON = HUGEFUN

      integer, parameter :: INT32 = selected_int_kind(8)
      ! INT32 is the kind for 32-bit integer, i.e., integer*4.
      ! SELECTED_INT_KIND(p) returns an INTEGER that equals the
      ! processor dependent kind type parameter of the integer type
      ! accommodating all values n with -10^p < n < 10^p. Therefore,
      ! SELECTED_INT_KIND(p) should be the kind of integer*4 for p = 7,
      ! 8, 9 (also for p = 5 and 6 unless integer*3 is supportred).

      end module consts 
