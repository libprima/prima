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

#include "ppf.h"

      module consts_mod 

#if __FORTRAN_STANDARD__ >= 2008
      use iso_fortran_env, only :INT16,INT32,INT64,REAL32,REAL64,REAL128
#endif

      implicit none
      private
      public :: DEBUG_MODE
      public :: IK, INT16, INT32, INT64, IK_DEFAULT
      public :: RP, DP, SP, QP, REAL32, REAL64, REAL128, RP_DEFAULT 
      public :: ZERO, ONE, TWO, HALF, QUART, TEN, TENTH, PI
      public :: EPS, HUGENUM, ALMOST_INFINITY, HUGEFUN, HUGECON 
      public :: SRNLEN

#if __DEBUG_MODE__ == 1
      logical, parameter :: DEBUG_MODE = .true.
#else
      logical, parameter :: DEBUG_MODE = .false.
#endif

#if __FORTRAN_STANDARD__ < 2008
      integer, parameter :: INT16 = selected_int_kind(4)
      integer, parameter :: INT32 = selected_int_kind(8)
      integer, parameter :: INT64 = selected_int_kind(16)
      integer, parameter :: REAL32 = selected_real_kind(6, 37)
      integer, parameter :: REAL64 = selected_real_kind(15, 307)
      integer, parameter :: REAL128 = selected_real_kind(33, 4931)
#endif      
      integer, parameter :: IK_DEFAULT = kind(0)  ! Default integer kind
      integer, parameter :: RP_DEFAULT = kind(0.0)  ! Default real kind
      integer, parameter :: SP = REAL32  ! Kind for single precision
      integer, parameter :: DP = REAL64  ! Kind for double precision
      integer, parameter :: QP = REAL128  ! Kind for quadruple precision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the integer kind to be used in the Fortran code.
#if __INTEGER_KIND__ == 16
      integer, parameter :: IK = INT16
#elif __INTEGER_KIND__ == 32 
      integer, parameter :: IK = INT32
#elif __INTEGER_KIND__ == 64 
      integer, parameter :: IK = INT64
#else
      integer, parameter :: IK = IK_DEFAULT
#endif
! Define the real kind to be used in the Fortran code.
#if __REAL_PRECISION__ == 32
      integer, parameter :: RP = REAL32  
#elif __REAL_PRECISION__ == 64 
      integer, parameter :: RP = REAL64  
#elif __REAL_PRECISION__ == 128 
      integer, parameter :: RP = REAL128
#else 
      integer, parameter :: RP = RP_DEFAULT
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      real(RP), parameter :: ZERO = 0.0_RP
      real(RP), parameter :: ONE = 1.0_RP
      real(RP), parameter :: TWO = 2.0_RP
      real(RP), parameter :: HALF = 0.5_RP
      real(RP), parameter :: QUART = 0.25_RP
      real(RP), parameter :: TEN = 10.0_RP
      real(RP), parameter :: TENTH = 0.1_RP
      real(RP), parameter ::PI=3.141592653589793238462643383279502884_RP
      ! We may set PI to acos(-1.0_RP), but F95 will complain about 
      ! "Elemental function as initialization expression with
      ! non-integer/non-character arguments".

      real(RP), parameter :: EPS = epsilon(ZERO)
      real(RP), parameter :: HUGENUM = huge(ZERO)
      real(RP), parameter :: ALMOST_INFINITY = HALF*HUGENUM
      
      integer, parameter :: MAXE = maxexponent(ZERO)
      real(RP), parameter :: HUGEFUN = TWO**min(100, MAXE/2)
      real(RP), parameter :: HUGECON = HUGEFUN

      ! The maximal length of subroutine/function names
      integer, parameter :: SRNLEN = 30


      end module consts_mod 
