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
 * FORTRAN_STANDARD_            which Fortran standard to follow: 2008, 2018, 2023
 * RELEASED_                    released or not: 1, 0
 * DEBUGGING_                   debug or not: 0, 1
 * INTEGER_KIND_                the integer kind to be used: 0, 32, 64, 16
 * REAL_PRECISION_              the real precision to be used: 64, 32, 128, 0
 * QP_AVAILABLE_                quad precision available or not: 0, 1
 * MAX_HIST_MEM_MB_             maximal MB memory for computation history: 300
 * AGGRESSIVE_OPTIONS_          compile the code with aggressive options: 0, 1
 *
 * N.B.:
 *
 * 0. USE THE DEFAULT IF UNSURE.
 *
 * 1. The macros can be modified by the -D option of the compilers. For example,
 *    -DDEBUGGING_=1 will set DEBUGGING_ to 1.
 *
 * 2. All the macros defined here have a trailing underscore. We avoid the following
 * patterns, because C and C++ reserve them for the implementation of the languages.
 * - Begins with two underscores
 * - Begins with underscore and uppercase letter
 * - Begins with underscore and something else
 * - Contains two consecutive underscores
 * See https://devblogs.microsoft.com/oldnewthing/20230109-00/?p=107685 for details.
 *
 * 3. Setting MAX_HIST_MEM_MB_ to a big value may lead failures due to large arrays.
 *
 * 4. If you change these macros, make sure that your compiler is supportive when
 * changing INTEGER_KIND_ and REAL_PRECISION_.
 *
 * 5. Why not define these macros as parameters in the Fortran code, e.g.,
 *
 * logical, parameter :: DEBUGGING_ == .false. ?
 *
 * Such a definition will work for DEBUGGING_, but not for the macros that
 * depend on the compiler. In addition, we can change the value of the macros
 * by the -D option of the compilers, which is impossible if we code them as
 * parameters in the Fortran code.
 *
 */
/******************************************************************************/


/******************************************************************************/
#if !defined PPF_H_  /* include guard to avoid double inclusion */
#define PPF_H_
/******************************************************************************/


/******************************************************************************/
/* Which Fortran standard to follow?
 * N.B.: The value of FORTRAN_STANDARD_ is NOT used in the code. We define it
 * only for the purpose of a record. */
#if !defined FORTRAN_STANDARD_
#define FORTRAN_STANDARD_ 2008 /* Will be default to 2018 later (in 2028?).   */
#endif
/******************************************************************************/


/******************************************************************************/
/* Is this a released version? Should be 1 except for the developers. */
#if !defined RELEASED_
#define RELEASED_ 1
#endif
/******************************************************************************/


/******************************************************************************/
/* Are we debugging?
 * RELEASED_ == 1 and DEBUGGING_ == 1 do not conflict. Users may debug. */
#if !defined DEBUGGING_
#define DEBUGGING_ 0
#endif
/******************************************************************************/


/******************************************************************************/
/* Which integer kind to use?
 * 0 = default INTEGER, 16 = INTEGER*2, 32 = INTEGER*4, 64 = INTEGER*8.
 * Make sure that your compiler supports the selected kind. */
#if !defined INTEGER_KIND_
#define INTEGER_KIND_ 0
#endif
/* Fortran standards guarantee that 0 is supported, but not the others. */
/******************************************************************************/


/******************************************************************************/
/* Which real kind to use?
 * 0 = default REAL (SINGLE PRECISION), 32 = REAL*4, 64 = REAL*8, 128 = REAL*16.
 * Make sure that your compiler supports the selected kind.  Note the following:
 * 1. The default REAL (i.e., 0) is the single-precision REAL.
 * 2. Fortran standards guarantee that 0, 32, and 64 are supported, but not 128.
 * 3. If you set REAL_PRECISION_ to 128, you must set QP_AVAILABLE_ to 1. */
#if !defined REAL_PRECISION_
#define REAL_PRECISION_ 64
#endif

/* Is quad precision available on this platform (compiler, hardware ...)? */
/* Note:
 * 1. Not all platforms support REAL128. For example, pgfortran 19 does not.
 * 2. It is not guaranteed that REAL128 has a wider range than REAL64. For
 *    example, REAL128 of nagfor 7.0 has a range of 291, while REAL64
 *    has a range of 307.
 * 3. It is rarely a good idea to use REAL128 as the working precision,
 *    which is probably inefficient and unnecessary.
 * 4. Set QP_AVAILABLE_ to 1 and REAL_PRECISION_ to 128 if REAL128 is available
 *    and you REALLY intend to use it. DO NOT DO THIS UNLESS REALLY SURE. */
#if !defined QP_AVAILABLE_
#define QP_AVAILABLE_ 0
#endif

/* Revise REAL_PRECISION_ according to QP_AVAILABLE_ . */
#if QP_AVAILABLE_ != 1 && REAL_PRECISION_ > 64
#undef REAL_PRECISION_
#define REAL_PRECISION_ 64
#endif
/******************************************************************************/


/******************************************************************************/
/* The maximal memory for recording the computation history (MB).
 * The maximal supported value is 2000, as 2000 M = 2*10^9 = maximum of INT32.
 * N.B.: A big value (even < 2000) may lead to SEGFAULTs due to large arrays. */
#if !defined MAX_HIST_MEM_MB_
#define MAX_HIST_MEM_MB_ 300  /* 1MB > 10^5*REAL64. 100 is sometimes too small.*/
#endif
/******************************************************************************/


/******************************************************************************/
/* Will we compile the code with aggressive options (e.g., -Ofast for gfortran)?
 * Some debugging will be disabled if yes (1). Note:
 * 1. It is OK to set AGGRESSIVE_OPTIONS_ = 0 and DEBUGGING_ = 1 simultaneously.
 * 2. When compiled with aggressive options, the code may behave unexpectedly.*/
#if !defined AGGRESSIVE_OPTIONS_
#define AGGRESSIVE_OPTIONS_ 0
#endif
/******************************************************************************/


/******************************************************************************/
#endif  /* include guard to avoid double inclusion */
/******************************************************************************/
