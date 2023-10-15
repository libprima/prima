#include "ppf.h"

module consts_mod
!--------------------------------------------------------------------------------------------------!
! This is a module defining some constants.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Sunday, October 15, 2023 PM12:03:14
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
! 3. Fortran standard (as of F2018) specifies the following for types INTEGER and REAL.
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
!    - Default integer, default real, and default logical all occupy one storage unit. Double
!    precision and complex occupy two storage units and double complex requires four storage units.
!
!    In other words, the standard only imposes that the following three types should be supported:
!    - INTEGER(KIND(0)), i.e., default integer,
!    - REAL(KIND(0.0)), i.e., default real (single-precision real),
!    - REAL(KIND(0.0D0)), i.e., double-precision real.
!
!    Moreover, the following should be noted.
!
!    - Other types of INTEGER/REAL may not be available on all platforms (e.g., nvfortran 23.3 and
!    flang 15.0.3 do not support REAL128).
!    - The standard does not specify the range of the default integer. However, if the default real
!    occupies 32 bits, which is normally the case, then the default integer occupies also 32 bits,
!    and hence the range is probably [-2^32, 2^31-1], approximately [-2*10^9, 2*10^9].
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

! Integer and real kinds. Unsupported kinds are negative.
use, intrinsic :: iso_fortran_env, only : INT16, INT32, INT64, REAL32, REAL64, REAL128

! Standard IO units
use, intrinsic :: iso_fortran_env, only : STDIN => INPUT_UNIT, STDOUT => OUTPUT_UNIT, STDERR => ERROR_UNIT

implicit none

private
public :: DEBUGGING
public :: IK, IK_DFT, INT16, INT32, INT64
public :: RP, RP_DFT, DP, SP, QP
public :: ZERO, ONE, TWO, HALF, QUART, TEN, TENTH, PI
public :: REALMIN, EPS, TINYCV, REALMAX, FUNCMAX, CONSTRMAX, BOUNDMAX
public :: SYMTOL_DFT
public :: STDIN, STDOUT, STDERR
public :: RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, CTOL_DFT, CWEIGHT_DFT
public :: ETA1_DFT, ETA2_DFT, GAMMA1_DFT, GAMMA2_DFT
public :: MAXFUN_DIM_DFT, MAXHISTMEM, MIN_MAXFILT, MAXFILT_DFT, IPRINT_DFT


logical, parameter :: DEBUGGING = (PRIMA_DEBUGGING == 1)  ! Whether we are in debugging mode
integer, parameter :: IK_DFT = kind(0)  ! Default integer kind
integer, parameter :: RP_DFT = kind(0.0)  ! Default real kind
integer, parameter :: SP = REAL32  ! Kind for single precision
integer, parameter :: DP = REAL64  ! Kind for double precision
integer, parameter :: QP = REAL128  ! Kind for quadruple precision

! Define the integer kind to be used in the Fortran code.
#if PRIMA_INTEGER_KIND  == 0
integer, parameter :: IK = IK_DFT
#elif PRIMA_INTEGER_KIND == 16
integer, parameter :: IK = INT16
#elif PRIMA_INTEGER_KIND == 32
integer, parameter :: IK = INT32
#elif PRIMA_INTEGER_KIND == 64
integer, parameter :: IK = INT64
#else
integer, parameter :: IK = IK_DFT
#endif

! Define the real kind to be used in the Fortran code.
#if PRIMA_REAL_PRECISION == 0
integer, parameter :: RP = RP_DFT
#elif PRIMA_REAL_PRECISION == 32
integer, parameter :: RP = REAL32
#elif PRIMA_REAL_PRECISION == 64
integer, parameter :: RP = REAL64
#elif PRIMA_REAL_PRECISION == 128
integer, parameter :: RP = REAL128
#else
integer, parameter :: RP = REAL64  ! Double precision
#endif

! Define some frequently used numbers.
real(RP), parameter :: ZERO = 0.0_RP
real(RP), parameter :: ONE = 1.0_RP
real(RP), parameter :: TWO = 2.0_RP
real(RP), parameter :: HALF = 0.5_RP
real(RP), parameter :: QUART = 0.25_RP
real(RP), parameter :: TEN = 10.0_RP
real(RP), parameter :: TENTH = 0.1_RP
real(RP), parameter :: PI = 3.1415926535897932384626433832795028841971693993751058209749445923078_RP

! EPS is the machine epsilon, namely the smallest floating-point number such that 1.0 + EPS > 1.0.
real(RP), parameter :: EPS = epsilon(ZERO)
! REALMIN is the smallest positive normalized floating-point number, which is 2^(-1022) ~ 2.225E-308
! for IEEE double precision. Taking double precision as an example, REALMIN in other languages:
! MATLAB: realmin or realmin('double')
! Python: numpy.finfo(numpy.float64).tiny
! Julia: realmin(Float64)
! R: double.xmin
real(RP), parameter :: REALMIN = tiny(ZERO)
! REALMAX is the largest positive floating-point number, which is 2^1023 * (2 - EPS) ~ 1.797E308
! for IEEE double precision. Taking double precision as an example, REALMAX in other languages:
! MATLAB: realmax or realmax('double')
! Python: numpy.finfo(numpy.float64).max
! Julia: realmax(Float64)
! R: double.xmax
real(RP), parameter :: REALMAX = huge(ZERO)

integer, parameter :: MINE = minexponent(ZERO)
integer, parameter :: MAXE = maxexponent(ZERO)

! TINYCV is used in LINCOA. Powell set TINYCV = 1.0D-60. What about setting TINYCV = REALMIN?
real(RP), parameter :: TINYCV = real(radix(ZERO), RP)**max(-200, MINE)  ! Normally, RADIX = 2.
! FUNCMAX is used in the moderated extreme barrier. All function values are projected to the
! interval [-FUNCMAX, FUNCMAX] before passing to the solvers, and NaN is replaced with FUNCMAX.
! CONSTRMAX plays a similar role for constraints.
real(RP), parameter :: FUNCMAX = real(radix(ZERO), RP)**min(100, MAXE / 2)  ! Normally, RADIX = 2.
real(RP), parameter :: CONSTRMAX = FUNCMAX
! Any bound with an absolute value at least BOUNDMAX is considered as no bound.
real(RP), parameter :: BOUNDMAX = QUART * REALMAX

! SYMTOL_DFT is the default tolerance for testing symmetry of matrices. It can be set to 0 if the
! IEEE Standard for Floating-Point Arithmetic (IEEE 754) is respected, particularly if addition and
! multiplication are commutative. However, as of 20220408, NAG nagfor does not ensure commutativity
! for REAL128. Indeed, Fortran standards do not enforce IEEE 754, so compilers are not guaranteed to
! respect it. Hence we set SYMTOL_DFT to a nonzero number when PRIMA_RELEASED is 1, although we do not
! intend to test symmetry in production. We set SYMTOL_DFT in the same way when PRIMA_DEBUGGING is 0.
! Update 20221226: When gfortran 12 is invoked with aggressive optimization options, it is buggy
! with ALL() and ANY(). We set SYMTOL_DFT to REALMAX to signify this case and disable the check.
! Update 20221229: ifx 2023.0.0 20221201 cannot ensure symmetry even up to 100*EPS if invoked
! with aggressive optimization options and if the floating-point numbers are in single precision.
! Update 20230307: ifx 2023.0.0 20221201 cannot ensure symmetry even up to 10*EPS if invoked with
! -O3 and if the floating-point numbers are in single precision.
! Update 20231002: HUAWEI BiSheng Compiler 2.1.0.B010 (flang) cannot ensure symmetry even up to
! 1.0E2*EPS if invoked with -Ofast and if the floating-point numbers are in single precision.
! This same is observed for arm-linux-compiler-22.1 on Kunpeng.
#if (defined __GFORTRAN__ || defined __INTEL_COMPILER && PRIMA_REAL_PRECISION < 64) && PRIMA_AGGRESSIVE_OPTIONS == 1
real(RP), parameter :: SYMTOL_DFT = REALMAX
#elif (defined __FLANG && PRIMA_REAL_PRECISION < 64) && PRIMA_AGGRESSIVE_OPTIONS == 1
real(RP), parameter :: SYMTOL_DFT = max(1.0E3 * EPS, 1.0E-10_RP)
#elif (defined __INTEL_COMPILER && PRIMA_REAL_PRECISION < 64)
real(RP), parameter :: SYMTOL_DFT = max(5.0E1 * EPS, 1.0E-10_RP)
#elif (defined __NAG_COMPILER_RELEASE && PRIMA_REAL_PRECISION > 64) || (PRIMA_RELEASED == 1) || (PRIMA_DEBUGGING == 0)
real(RP), parameter :: SYMTOL_DFT = max(1.0E1 * EPS, 1.0E-10_RP)
#else
real(RP), parameter :: SYMTOL_DFT = ZERO
#endif

! Some default values
real(RP), parameter :: RHOBEG_DFT = ONE
real(RP), parameter :: RHOEND_DFT = 1.0E-6_RP
real(RP), parameter :: FTARGET_DFT = -REALMAX
real(RP), parameter :: CTOL_DFT = sqrt(EPS)
real(RP), parameter :: CWEIGHT_DFT = 1.0E8_RP
real(RP), parameter :: ETA1_DFT = TENTH
real(RP), parameter :: ETA2_DFT = 0.7_RP
real(RP), parameter :: GAMMA1_DFT = HALF
real(RP), parameter :: GAMMA2_DFT = TWO
integer(IK), parameter :: IPRINT_DFT = 0
integer(IK), parameter :: MAXFUN_DIM_DFT = 500

! Maximal amount of memory (Byte) allowed for XHIST, FHIST, CONHIST, CHIST, and the filters.
integer, parameter :: MHM = PRIMA_MAX_HIST_MEM_MB * 10**6
! Make sure that MAXHISTMEM does not exceed HUGE(0) to avoid overflow and memory errors.
integer, parameter :: MAXHISTMEM = min(MHM, huge(0))

! Maximal length of the filter used in constrained solvers.
integer(IK), parameter :: MIN_MAXFILT = 200  ! Should be positive; < 200 is not recommended.
integer(IK), parameter :: MAXFILT_DFT = 10_IK * MIN_MAXFILT


end module consts_mod
