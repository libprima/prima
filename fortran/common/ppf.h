/******************************************************************************/
/*
 * Coded by Zaikun ZHANG (www.zhangzk.net) in July 2020.
 *
 * Last Modified: Tue 11 April 2022 04:53:00 PM HKT
 */
/******************************************************************************/
/*
 * ppf.h defines the following preprocessing macros (the first value is default).
 *
 * PRIMA_FORTRAN_STANDARD   which Fortran standard to follow: 2008, 2018, 2023
 * PRIMA_RELEASED           released or not: 1, 0
 * PRIMA_DEBUGGING          debug or not: 0, 1
 * PRIMA_INTEGER_KIND       the integer kind to be used: 0, 32, 64, 16
 * PRIMA_REAL_PRECISION     the real precision to be used: 64, 32, 128, 0
 * PRIMA_QP_AVAILABLE       quad precision available or not: 0, 1
 * PRIMA_MAX_HIST_MEM_MB    maximal MB memory for computation history: 300
 * PRIMA_AGGRESSIVE_OPTIONS compile the code with aggressive options: 0, 1
 *
 * N.B.:
 *
 * 0. USE THE DEFAULT IF UNSURE.
 *
 * 1. The macros can be modified by the -D option of the compilers. For example,
 *    -DPRIMA_DEBUGGING=1 will set PRIMA_DEBUGGING to 1.
 *
 * 2. All the macros defined here starts with "PRIMA_". We avoid the following
 *    patterns, as C and C++ reserve them for the implementation of the languages.
 *    - Begins with two underscores
 *    - Begins with underscore and uppercase letter
 *    - Begins with underscore and something else
 *    - Contains two consecutive underscores
 *   See https://devblogs.microsoft.com/oldnewthing/20230109-00/?p=107685.
 *
 * 3. Setting PRIMA_MAX_HIST_MEM_MB to a big value may lead to segfaults due to
 *    large arrays.
 *
 * 4. If you change these macros, make sure that your compiler is supportive
 *    when changing PRIMA_INTEGER_KIND and PRIMA_REAL_PRECISION.
 *
 * 5. Why not define these macros as parameters in the Fortran code, e.g.,
 *
 *    logical, parameter :: PRIMA_DEBUGGING == .false. ?
 *
 *    Such a definition will work for PRIMA_DEBUGGING, but not for the macros
 *    that depend on the compiler. In addition, we can change the value of the
 *    macros by the -D option of the compilers, which is impossible if we code
 *    them as parameters in the Fortran code.
 *
 */
/******************************************************************************/


/******************************************************************************/
#if !defined PRIMA_PPF_H  /* include guard to avoid double inclusion */
#define PRIMA_PPF_H
/******************************************************************************/


/******************************************************************************/
/* Which Fortran standard to follow?
 * N.B.: 1. The value of PRIMA_FORTRAN_STANDARD is NOT used in the code. We
 * define it only for the purpose of a record.
 * 2. With gfortran, due to `error stop` and `backtrace`, we must either compile
 * with no `-std` or use `-std=f20xy -fall-intrinsics` with xy >= 18. */
#if !defined PRIMA_FORTRAN_STANDARD
#define PRIMA_FORTRAN_STANDARD 2008  /* Default to 2018 later (in 2025?).     */
#endif
/******************************************************************************/


/******************************************************************************/
/* Is this a released version? Should be 1 except for the developers. */
#if !defined PRIMA_RELEASED
#define PRIMA_RELEASED 1
#endif
/******************************************************************************/


/******************************************************************************/
/* Are we debugging?
 * PRIMA_RELEASED == 1 and PRIMA_DEBUGGING == 1 do not conflict. User may debug.*/
#if !defined PRIMA_DEBUGGING
#define PRIMA_DEBUGGING 0
#endif
/******************************************************************************/


/******************************************************************************/
/* Which integer kind to use?
 * 0 = default INTEGER, 16 = INTEGER*2, 32 = INTEGER*4, 64 = INTEGER*8.
 * Make sure that your compiler supports the selected kind. */
#if !defined PRIMA_INTEGER_KIND
#define PRIMA_INTEGER_KIND 0
#endif
/* Fortran standards guarantee that 0 is supported, but not the others. */
/******************************************************************************/


/******************************************************************************/
/* Which real kind to use?
 * 0 = default REAL (SINGLE PRECISION), 32 = REAL*4, 64 = REAL*8, 128 = REAL*16.
 * Make sure that your compiler supports the selected kind.  Note the following:
 * 1. The default REAL (i.e., 0) is the single-precision REAL.
 * 2. Fortran standards guarantee that 0, 32, and 64 are supported, but not 128.
 * 3. If you set PRIMA_REAL_PRECISION to 128, then set PRIMA_QP_AVAILABLE to 1.*/
#if !defined PRIMA_REAL_PRECISION
#define PRIMA_REAL_PRECISION 64
#endif

/* Is quad precision available on this platform (compiler, hardware ...)? */
/* Note:
 * 1. Not all platforms support REAL128. For example, pgfortran 19 does not.
 * 2. It is not guaranteed that REAL128 has a wider range than REAL64. For
 *    example, REAL128 of nagfor 7.0 has a range of 291, while REAL64
 *    has a range of 307.
 * 3. It is rarely a good idea to use REAL128 as the working precision,
 *    which is probably inefficient and unnecessary.
 * 4. Set PRIMA_QP_AVAILABLE to 1 and PRIMA_REAL_PRECISION to 128 if REAL128
 *    is available and you REALLY intend to use it. DO NOT DO IT IF NOT SURE. */
#if !defined PRIMA_QP_AVAILABLE
#define PRIMA_QP_AVAILABLE 0
#endif

/* Revise PRIMA_REAL_PRECISION according to PRIMA_QP_AVAILABLE . */
#if PRIMA_QP_AVAILABLE != 1 && PRIMA_REAL_PRECISION > 64
#undef PRIMA_REAL_PRECISION
#define PRIMA_REAL_PRECISION 64
#endif
/******************************************************************************/


/******************************************************************************/
/* The maximal memory for recording the computation history (MB).
 * The maximal supported value is 2000, as 2000 M = 2*10^9 = maximum of INT32.
 * N.B.: A big value (even < 2000) may lead to SEGFAULTs due to large arrays. */
#if !defined PRIMA_MAX_HIST_MEM_MB
#define PRIMA_MAX_HIST_MEM_MB 300  /* 1MB > 10^5*REAL64. 100 can be too small.*/
#endif
/******************************************************************************/


/******************************************************************************/
/* Will we compile the code with aggressive options (e.g., -Ofast for gfortran)?
 * Some debugging will be disabled if yes (1). Note:
 * 1. It is OK to set PRIMA_AGGRESSIVE_OPTIONS = 0 and PRIMA_DEBUGGING = 1
 *    simultaneously.
 * 2. When compiled with aggressive options, the code may behave unexpectedly.*/
#if !defined PRIMA_AGGRESSIVE_OPTIONS
#define PRIMA_AGGRESSIVE_OPTIONS 0
#endif
/******************************************************************************/


/******************************************************************************/
#endif  /* include guard to avoid double inclusion */
/******************************************************************************/
