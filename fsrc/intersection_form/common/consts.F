!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of consts.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 11-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! CONSTS_MOD is a module defining some constants.
!
! Coded by Zaikun Zhang in July 2020.

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
      use intrinsic :: iso_fortran_env, only : INT16, INT32, INT64, REAL&
     &32, REAL64, REAL128
#endif

      implicit none
      private
      public :: DEBUGGING
      public :: IK, INT16, INT32, INT64, IK_DFT
      public :: RP, DP, SP, QP, REAL32, REAL64, REAL128, RP_DFT
      public :: ZERO, ONE, TWO, HALF, QUART, TEN, TENTH, PI
      public :: EPS, HUGENUM, ALMOST_INFINITY, HUGEFUN, HUGECON
      public :: SRNLEN, MSSGLEN
      public :: OUTUNIT
      public :: RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, IPRINT_DFT
      public :: MAXFUN_DIM_DFT, MAXIMAL_HIST

#if __DEBUGGING__ == 1
      logical, parameter :: DEBUGGING = .true.
#else
      logical, parameter :: DEBUGGING = .false.
#endif

#if __FORTRAN_STANDARD__ < 2008
! With gfortran: selected_real_kind(k) will return INT16 with k = 3 and
! 4, return INT32 with k = 5--9, and return INT64 with k = 10--18.
      integer, parameter :: INT16 = selected_int_kind(4)
      integer, parameter :: INT32 = selected_int_kind(7)
      integer, parameter :: INT64 = selected_int_kind(14)
      integer, parameter :: REAL32 = kind(0.0)
      integer, parameter :: REAL64 = kind(0.0D0)
      integer, parameter :: REAL128 = selected_real_kind(p = 30)
#endif
      integer, parameter :: IK_DFT = kind(0) ! Default integer kind
      integer, parameter :: RP_DFT = kind(0.0) ! Default real kind
      integer, parameter :: SP = REAL32 ! Kind for single precision
      integer, parameter :: DP = REAL64 ! Kind for double precision
      integer, parameter :: QP = REAL128 ! Kind for quadruple precision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the integer kind to be used in the Fortran code.
#if __INTEGER_KIND__ == 16
      integer, parameter :: IK = INT16
#elif __INTEGER_KIND__ == 32
      integer, parameter :: IK = INT32
#elif __INTEGER_KIND__ == 64
      integer, parameter :: IK = INT64
#else
      integer, parameter :: IK = IK_DFT
#endif
! Define the real kind to be used in the Fortran code.
#if __REAL_PRECISION__ == 32
      integer, parameter :: RP = REAL32
#elif __REAL_PRECISION__ == 64
      integer, parameter :: RP = REAL64
#elif __REAL_PRECISION__ == 128
      integer, parameter :: RP = REAL128
#else
      integer, parameter :: RP = REAL64 ! double precision
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(RP), parameter :: ZERO = 0.0_RP
      real(RP), parameter :: ONE = 1.0_RP
      real(RP), parameter :: TWO = 2.0_RP
      real(RP), parameter :: HALF = 0.5_RP
      real(RP), parameter :: QUART = 0.25_RP
      real(RP), parameter :: TEN = 10.0_RP
      real(RP), parameter :: TENTH = 0.1_RP
      real(RP), parameter ::PI = 3.141592653589793238462643383279502884_&
     &RP
! We may set PI to acos(-1.0_RP), but some compilers may complain about
! "Elemental function as initialization expression with non-integer or
! non-character arguments".

      real(RP), parameter :: EPS = epsilon(ZERO)
      real(RP), parameter :: HUGENUM = huge(ZERO)
      real(RP), parameter :: ALMOST_INFINITY = HALF*HUGENUM

      integer, parameter :: MAXE = maxexponent(ZERO)
      real(RP), parameter :: HUGEFUN = TWO**min(100, MAXE/2)
      real(RP), parameter :: HUGECON = HUGEFUN

! The maximal length of subroutine/function names
      integer, parameter :: SRNLEN = 30

! The maximal length of messages; used in output.F and mexapi.F
      integer, parameter :: MSSGLEN = 300

! Output unit, can be any integer between 9 and 99
      integer(IK), parameter :: OUTUNIT = 9

! Some default values
      real(RP), parameter :: RHOBEG_DFT = ONE
      real(RP), parameter :: RHOEND_DFT = 1.0E-6_RP
      real(RP), parameter :: FTARGET_DFT = -HUGENUM
      integer(IK), parameter :: IPRINT_DFT = 1
      integer(IK), parameter :: MAXFUN_DIM_DFT = 500

! Maximal amount of history that XHIST, FHIST, CHIST can return
#if __INTEGER_KIND__ <= 16
      integer(IK), parameter :: MAXIMAL_HIST = huge(0_IK)
#else
      integer(IK), parameter :: MAXIMAL_HIST = 100000
#endif


      end module consts_mod