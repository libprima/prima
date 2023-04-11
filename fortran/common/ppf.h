/******************************************************************************/
/*
 * Coded by Zaikun ZHANG (www.zhangzk.net) in July 2020.
 *
 * Last Modified: Tue 11 April 2022 04:53:00 PM HKT
 */
/******************************************************************************/
/*
 * ppf.h defines the following preprocessing flags (the first value is default).
 *
 * RELEASED_                    released or not: 1, 0
 * DEBUGGING_                   debug or not: 0, 1
 * FORTRAN_STANDARD_            which Fortran standard to follow: 2008, 2018, 2023
 * INTEGER_KIND_                the integer kind to be used: 0, 32, 64, 16
 * REAL_PRECISION_              the real precision to be used: 64, 32, 128, 0
 * QP_AVAILABLE_                quad precision available or not: 0, 1
 * MAX_HIST_MEM_MB_             maximal MB memory for computation history: 300
 * AGGRESSIVE_OPTIONS_          compile the code with aggressive options: 0, 1
 * USE_STORAGE_SIZE_            use the STORAGE_SIZE intrinsic or not: 0, 1
 * USE_ISO_FORTRAN_ENV_INTREAL_ use INT32 etc in ISO_FORTRAN_ENV or not: 0, 1
 *
 * N.B.:
 *
 * 0. USE THE DEFAULT IF UNSURE.
 *
 * 1. All the macros defined here have a trailing underscore. We avoid the following
 * patterns, because C and C++ reserve them for the implementation of the languages.
 * - Begins with two underscores
 * - Begins with underscore and uppercase letter
 * - Begins with underscore and something else
 * - Contains two consecutive underscores
 * See https://devblogs.microsoft.com/oldnewthing/20230109-00/?p=107685 for details.
 *
 * 2. Setting MAX_HIST_MEM_MB_ to a big value may lead failures due to large arrays.
 *
 * 3. If you change these flags, make sure that your compiler is supportive when
 * changing INTEGER_KIND_, REAL_PRECISION_, FORTRAN_STANDARD_, USE_STORAGE_SIZE_
 * (Fortran 2008), USE_ISO_FORTRAN_ENV_INTREAL_ (Fortran 2008).
 *
 * 4. Later, when Fortran 2008 is more widely supported by the compilers (e.g.,
 * in 2025?), we will default FORTRAN_STANDARD_ to 2008. In addition, we will
 * remove USE_STORAGE_SIZE_ and USE_ISO_FORTRAN_ENV_INTREAL_ and modify
 * the package so that everything behaves as if the flags are both 1.
 *
 * 5. Why not define these flags as parameters in the Fortran code, e.g.,
 *
 * logical, parameter :: DEBUGGING_ == .false. ?
 *
 * Such a definition will work for DEBUGGING_, but not for the flags that
 * depend on the compiler, for instance, USE_ISO_FORTRAN_ENV_INTREAL_.
 *
 */
/******************************************************************************/


/******************************************************************************/
#if !defined PPF_H_  /* include guard to avoid double inclusion */
#define PPF_H_
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
/* Which Fortran standard to follow? */
#if !defined FORTRAN_STANDARD_
#define FORTRAN_STANDARD_ 2008 /* Will be default to 2018 later (in 2028?).*/
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
 * 2. If you set REAL_PRECISION_ to 128, you must set QP_AVAILABLE_ to 1. */
#if !defined REAL_PRECISION_
#define REAL_PRECISION_ 64
#endif
/* Fortran standards guarantee that 0, 32, and 64 are supported, but not 128. */

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
#define MAX_HIST_MEM_MB_ 300  /* 1MB > 10^5*REAL64. 100 is sometimes too small. */
#endif
/******************************************************************************/


/******************************************************************************/
/* Will we compile the code with aggressive options (e.g., -Ofast for gfortran)?
 * Some debugging will be disabled if yes (1). Note:
 * 1. It is VALID to set AGGRESSIVE_OPTIONS_ to 0 and DEBUGGING_ to 1 at
 * the same time.
 * 2. When compiled with aggressive options, the code may behave unexpectedly */
#if !defined AGGRESSIVE_OPTIONS_
#define AGGRESSIVE_OPTIONS_ 0
#endif
/******************************************************************************/


/******************************************************************************/
/* Do we use the STORAGE_SIZE intrinsic? (Fortran 2008) */
/* We prefer STORAGE_SIZE to C_SIZEOF, because the former is intrinsic while the
 * later requires the intrinsic module ISO_C_BINDING.
 * The flag will be removed later (in 2025?). */
#if !defined USE_STORAGE_SIZE_
#define USE_STORAGE_SIZE_ 0
#endif

#if defined __GFORTRAN__
#if __GNUC__ >= 5 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#undef USE_STORAGE_SIZE_
#define USE_STORAGE_SIZE_ 1
#endif
#endif

#if defined __INTEL_COMPILER
#if __INTEL_COMPILER >= 1710
#undef USE_STORAGE_SIZE_
#define USE_STORAGE_SIZE_ 1
#endif
#endif

#if defined __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE >= 53
#undef USE_STORAGE_SIZE_
#define USE_STORAGE_SIZE_ 1
#endif
#endif

#if defined __PGI
#if __PGIC__ >= 15 && __PGIC_MINOR__ >= 4
#undef USE_STORAGE_SIZE_
#define USE_STORAGE_SIZE_ 1
#endif
#endif

#if defined __ibmxl__
#if __ibmxl_version__ >= 15 && __ibmxl_release__ >= 2
#undef USE_STORAGE_SIZE_
#define USE_STORAGE_SIZE_ 1
#endif
#endif

#if FORTRAN_STANDARD_ < 2008
#undef USE_STORAGE_SIZE_
#define USE_STORAGE_SIZE_ 0
#endif
/******************************************************************************/


/******************************************************************************/
/* Do we use INT16, INT32, INT64, REAL32, REAL64, REAL128 from ISO_FORTRAN_ENV?
 * (Fortran 2008). This flag will be removed later (in 2025?). */
#if !defined USE_ISO_FORTRAN_ENV_INTREAL_
#define USE_ISO_FORTRAN_ENV_INTREAL_ 0
#endif

#if defined __GFORTRAN__
#if __GNUC__ >= 5 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5)
#undef USE_ISO_FORTRAN_ENV_INTREAL_
#define USE_ISO_FORTRAN_ENV_INTREAL_ 1
#endif
#endif

#if defined __INTEL_COMPILER
#if __INTEL_COMPILER >= 1640
#undef USE_ISO_FORTRAN_ENV_INTREAL_
#define USE_ISO_FORTRAN_ENV_INTREAL_ 1
#endif
#endif

#if defined __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE >= 53
#undef USE_ISO_FORTRAN_ENV_INTREAL_
#define USE_ISO_FORTRAN_ENV_INTREAL_ 1
#endif
#endif

#if defined __PGI
#if __PGIC__ >= 14 && __PGIC_MINOR__ >= 1
#undef USE_ISO_FORTRAN_ENV_INTREAL_
#define USE_ISO_FORTRAN_ENV_INTREAL_ 1
#endif
#endif

#if defined __ibmxl__
#if __ibmxl_version__ >= 14 && __ibmxl_release__ >= 1
#undef USE_ISO_FORTRAN_ENV_INTREAL_
#define USE_ISO_FORTRAN_ENV_INTREAL_ 1
#endif
#endif

#if FORTRAN_STANDARD_ < 2008
#undef USE_ISO_FORTRAN_ENV_INTREAL_
#define USE_ISO_FORTRAN_ENV_INTREAL_ 0
#endif
/******************************************************************************/


/******************************************************************************/
#endif  /* include guard to avoid double inclusion */
/******************************************************************************/
