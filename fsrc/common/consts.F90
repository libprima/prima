#include "ppf.h"

module consts_mod
!--------------------------------------------------------------------------------------------------!
! This is a module defining some constants.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Sunday, December 19, 2021 AM11:30:27
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Remarks:
!
! 1. REAL*4, REAL*8, INTEGER*4, INTEGER*8 are not Fortran standard expressions. Do not use them!
!
! 2. Never use KIND with a literal value, e.g., REAL(KIND = 8), because Fortran standards never
! define what KIND = 8 means. There is NO guarantee that REAL(KIND = 8) will be legal, let alone
! being double precision.
!
! 3. Fortran standard (as of F2003) specifies the following for types INTEGER and REAL.
!
!    - A processor shall provide ONE OR MORE representation methods that define sets of values for
!    data of type integer; if the kind type parameter is not specified, the default kind value is
!    KIND(0) and the type specified is DEFAULT INTEGER.
!    - A processor shall provide TWO OR MORE approximation methods that define sets of values for
!    data of type real; if the type keyword REAL is specified and the kind type parameter is not
!    specified, the default kind value is KIND (0.0) and the type specified is DEFAULT REAL; If the
!    type keyword DOUBLE PRECISION is specified, the kind value is KIND (0.0D0) and the type
!    specified is DOUBLE PRECISION real; the decimal precision of the double precision real
!    approximation method shall be greater than that of the default real method.
!
!    In other words, the standard only imposes that the following three types should be supported:
!    - INTEGER(KIND(0)), i.e., default integer,
!    - REAL(KIND(0.0)), i.e., default real (single-precision real),
!    - REAL(KIND(0.0D0)), i.e., double-precision real.
!
!    Therefore, the following should be noted.
!
!    - Other types of INTEGER/REAL may not be available on all platforms (e.g., nvfortran 20 and
!    flang 7.1.0 do not support REAL128).
!    - The standard does not specify the range of the default integer.
!    - The standard does not specify what the range and precision of the default real or the
!    double-precision real, except that KIND(0.0D0) should have a greater precision than KIND(0.0)
!    --- no requirement about the range.
!
!    Consequently, the following should be observed in all Fortran code.
!
!    - DO NOT use any kind parameter other than IK, IK_DFT, RP, RP_DFT, SP, or DP, unless you are
!    sure that it is supported by your platform.
!    - DO NOT make any assumption on the range of INTEGER, REAL, or REAL(0.0D0) unless you are sure.
!    - Be cautious about OVERFLOW! In particular, for integers working as the lower/upper limit of
!    arrays, overflow can lead to Segmentation Faults!
!--------------------------------------------------------------------------------------------------!

#if __USE_ISO_FORTRAN_ENV_INTREAL__ == 1

#if __INTEGER_KIND__ == 16
use, intrinsic :: iso_fortran_env, only : INT16
#elif __INTEGER_KIND__ == 32
use, intrinsic :: iso_fortran_env, only : INT32
#elif __INTEGER_KIND__ == 64
use, intrinsic :: iso_fortran_env, only : INT64
#endif


use, intrinsic :: iso_fortran_env, only : REAL32, REAL64, REAL128
! The unsupported kind parameter will be negative.
#endif

implicit none
private
public :: DEBUGGING
public :: IK, IK_DFT
public :: RP, DP, SP, QP, RP_DFT
public :: ZERO, ONE, TWO, HALF, QUART, TEN, TENTH, PI
public :: REALMIN, EPS, HUGENUM, ALMOST_INFINITY, HUGEFUN, HUGECON
public :: MSSGLEN, FNAMELEN
public :: OUTUNIT
public :: RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, CTOL_DFT, IPRINT_DFT
public :: ETA1_DFT, ETA2_DFT, GAMMA1_DFT, GAMMA2_DFT
public :: MAXFUN_DIM_DFT, MAXMEMORY, MAXFILT_DFT


#if __DEBUGGING__ == 1
logical, parameter :: DEBUGGING = .true.
#else
logical, parameter :: DEBUGGING = .false.
#endif

#if __USE_ISO_FORTRAN_ENV_INTREAL__ != 1
! For gfortran, SELECTED_REAL_KIND(K) returns INT16 with K = 3--4, INT32 with k = 5--9, and INT64
! with K = 10--18. SELECTED_REAL_KIND returns a negative value for an unsupported kind.
#if __INTEGER_KIND__ == 16
integer, parameter :: INT16 = selected_int_kind(4)
#elif __INTEGER_KIND__ == 32
integer, parameter :: INT32 = selected_int_kind(7)
#elif __INTEGER_KIND__ == 64
integer, parameter :: INT64 = selected_int_kind(14)
#endif

integer, parameter :: REAL32 = kind(0.0)
integer, parameter :: REAL64 = kind(0.0D0)
integer, parameter :: REAL128 = selected_real_kind(p=30)

#endif
integer, parameter :: IK_DFT = kind(0)  ! Default integer kind
integer, parameter :: RP_DFT = kind(0.0)  ! Default real kind
integer, parameter :: SP = REAL32  ! Kind for single precision
integer, parameter :: DP = REAL64  ! Kind for double precision
integer, parameter :: QP = REAL128  ! Kind for quadruple precision

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the integer kind to be used in the Fortran code.
#if __INTEGER_KIND__  == 0
integer, parameter :: IK = IK_DFT
#elif __INTEGER_KIND__ == 16
integer, parameter :: IK = INT16
#elif __INTEGER_KIND__ == 32
integer, parameter :: IK = INT32
#elif __INTEGER_KIND__ == 64
integer, parameter :: IK = INT64
#else
integer, parameter :: IK = IK_DFT
#endif
! Define the real kind to be used in the Fortran code.
#if __REAL_PRECISION__ == 0
integer, parameter :: RP = RP_DFT
#elif __REAL_PRECISION__ == 32
integer, parameter :: RP = REAL32
#elif __REAL_PRECISION__ == 64
integer, parameter :: RP = REAL64
#elif __REAL_PRECISION__ == 128
integer, parameter :: RP = REAL128
#else
integer, parameter :: RP = REAL64  ! double precision
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(RP), parameter :: ZERO = 0.0_RP
real(RP), parameter :: ONE = 1.0_RP
real(RP), parameter :: TWO = 2.0_RP
real(RP), parameter :: HALF = 0.5_RP
real(RP), parameter :: QUART = 0.25_RP
real(RP), parameter :: TEN = 10.0_RP
real(RP), parameter :: TENTH = 0.1_RP
real(RP), parameter :: PI = 3.141592653589793238462643383279502884_RP
! We may set PI to acos(-1.0_RP), but some compilers may complain about `Elemental function as
! initialization expression with non-integer or non-character arguments`.

! REALMIN is the smallest positive normalized floating-point number, which is 2^(-1022), ~2.225E-308
! for IEEE double precision. Taking double precision as an example, REALMIN in other languages:
! MATLAB: realmin or realmin('double')
! Python: numpy.finfo(numpy.float64).tiny
! Julia: realmin(Float64)
real(RP), parameter :: REALMIN = tiny(ZERO)
real(RP), parameter :: EPS = epsilon(ZERO)  ! Machine epsilon
real(RP), parameter :: HUGENUM = huge(ZERO)
real(RP), parameter :: ALMOST_INFINITY = HALF * HUGENUM

integer, parameter :: MAXE = maxexponent(ZERO)
real(RP), parameter :: HUGEFUN = TWO**min(100, MAXE / 2)
real(RP), parameter :: HUGECON = HUGEFUN

! The maximal length of messages; used in output.f90 and fmexapi.F90
integer, parameter :: MSSGLEN = 1000

! The maximal length of output file names; used in output.f90
integer, parameter :: FNAMELEN = 1000

! Output unit, can be any integer between 9 and 99; used in output.f90
integer, parameter :: OUTUNIT = 9

! Some default values
real(RP), parameter :: RHOBEG_DFT = ONE
real(RP), parameter :: RHOEND_DFT = 1.0E-6_RP
real(RP), parameter :: FTARGET_DFT = -HUGENUM
real(RP), parameter :: CTOL_DFT = EPS
real(RP), parameter :: ETA1_DFT = TENTH
real(RP), parameter :: ETA2_DFT = 0.7_RP
real(RP), parameter :: GAMMA1_DFT = HALF
real(RP), parameter :: GAMMA2_DFT = TWO
integer(IK), parameter :: IPRINT_DFT = 0_IK
integer(IK), parameter :: MAXFUN_DIM_DFT = 500_IK

! Maximal amount of memory (Byte) allowed for XHIST, FHIST, CONHIST, CHIST, and the filters.
integer, parameter :: MXMMY = 21 * (10**8)   ! 21*10**8 = 2G.
! Make sure that MAXMEMORY does not exceed HUGE(0) to avoid overflow and memory errors.
integer, parameter :: MAXMEMORY = min(MXMMY, huge(0))

! Maximal length of the filter used in constrained solvers.
integer, parameter :: MAXFILT_DFT = 2000_IK

end module consts_mod
