/******************************************************************************/
/*
 * Coded by Zaikun ZHANG (www.zhangzk.net) in July 2020.
 *
 * Last Modified: Mon 24 May 2021 04:21:00 PM HKT
 */
/******************************************************************************/
/*
 * ppf.h defines the following preprocessing flags (the first value is default).
 *
 * __RELEASED__                released or not: 1, 0
 * __DEBUGGING__               debug or not: 0, 1
 * __FORTRAN_STANDARD__        which Fortran standard to follow: 2008, 2018, 2023
 * __INTEGER_KIND__            the integer kind to be used: 0, 32, 64, 16
 * __REAL_PRECISION__          the real precision to be used: 64, 32, 128, 0
 * __MAXHISTMEM__         maximal MB memory for computation history: 100
 * __AGRESSIVE_OPTIONS__       compile the code with aggressive options: 0, 1
 * __USE_STORAGE_SIZE__        use the STORAGE_SIZE intrinsic or not: 0, 1
 * __USE_ISO_FORTRAN_ENV_INTREAL__ use INT32 etc in ISO_FORTRAN_ENV or not: 0, 1
 *
 * N.B.:
 *
 * 0. USE THE DEFAULT IF UNSURE.
 *
 * 1. Setting __MAXHISTMEM__ to a big value may lead failures due to large arrays.
 *
 * 2. If you change these flags, make sure that your compiler is supportive
 * when changing __INTEGER_KIND__, __REAL_PRECISION__, __FORTRAN_STANDARD__,
 * __USE_STORAGE_SIZE__ (Fortran 2008),
 * __USE_ISO_FORTRAN_ENV_INTREAL__ (Fortran 2008).
 *
 * 3. Later, when Fortran 2008 is more widely supported by the compilers (e.g.,
 * in 2025?), we will default __FORTRAN_STANDARD__ to 2008. In addition, we will
 * remove __USE_STORAGE_SIZE__ and __USE_ISO_FORTRAN_ENV_INTREAL__ and modify
 * the package so that everything behaves as if the flags are both 1.
 *
 * 4. Why not define these flags as parameters in the Fortran code, e.g.,
 *
 * logical, parameter :: __DEBUGGING__ == .false. ?
 *
 * Such a definition will work for __DEBUGGING__, but not for the flags that
 * depend on the compiler, for instance, __USE_ISO_FORTRAN_ENV_INTREAL__.
 *
 */
/******************************************************************************/


/******************************************************************************/
/* Is this a released version? Should be 1 except for the developers. */
#if defined __RELEASED__
#undef __RELEASED__
#endif
#define __RELEASED__ 1
/******************************************************************************/


/******************************************************************************/
/* Are we debugging?
 * __RELEASED__ == 1 and __DEBUGGING__ == 1 do not conflict. Users may debug. */
#if defined __DEBUGGING__
#undef __DEBUGGING__
#endif
#define __DEBUGGING__ 0
/******************************************************************************/


/******************************************************************************/
/* Which Fortran standard to follow? */
#if defined __FORTRAN_STANDARD__
#undef __FORTRAN_STANDARD__
#endif
#define __FORTRAN_STANDARD__ 2008 /* Will be default to 2018 later (in 2028?).*/
/******************************************************************************/


/******************************************************************************/
/* Which integer kind to use?
 * 0 = default INTEGER, 16 = INTEGER*2, 32 = INTEGER*4, 64 = INTEGER*8.
 * Make sure that your compiler supports the selected kind. */
#if defined __INTEGER_KIND__
#undef __INTEGER_KIND__
#endif
#define __INTEGER_KIND__ 0
/* Fortran standards guarantee that 0 is supported, but not the others. */
/******************************************************************************/


/******************************************************************************/
/* Which real kind to use?
 * 0 = default REAL (SINGLE PRECISION), 32 = REAL*4, 64 = REAL*8, 128 = REAL*16.
 * Make sure that your compiler supports the selected kind.  Note the following:
 * 1. The default REAL (i.e., 0) is the single-precision REAL.
 * 2. If you set __REAL_PRECISION__ to 128, you must set __QP_AVAILABLE__ to 1. */
#if defined __REAL_PRECISION__
#undef __REAL_PRECISION__
#endif
#define __REAL_PRECISION__ 64
/* Fortran standards guarantee that 0, 32, and 64 are supported, but not 128. */

/* Is quad precision available on this platform (compiler, hardware ...)? */
/* Note:
 * 1. Not all platforms support REAL128. For example, pgfortran 19 does not.
 * 2. It is not guaranteed that REAL128 has a wider range than REAL64. For
 *    example, REAL128 of nagfor 7.0 has a range of 291, while REAL64
 *    has a range of 307.
 * 3. It is rarely a good idea to use REAL128 as the working precision,
 *    which is probably inefficient and unnecessary. */
#if defined __QP_AVAILABLE__
#undef __QP_AVAILABLE__
#endif
/* Change the following line to set __QP_AVAILABLE__ to 1 if REAL128 is
 * available and if you intend to use it. DO NOT CHANGE IT UNLESS REALLY SURE. */
#define __QP_AVAILABLE__ 0

/* Revise __REAL_PRECISION__ according to __QP_AVAILABLE__ . */
#if __QP_AVAILABLE__ != 1 && __REAL_PRECISION__ > 64
#undef __REAL_PRECISION__
#define __REAL_PRECISION__ 64
#endif
/******************************************************************************/


/******************************************************************************/
/* The maximal memory for recording the computation history (MB).
 * The maximal supported value is 2000, as 2000 M = 2*10^9 = maximum of INT32.
 * N.B.: A big value (even < 2000) may lead to SEGFAULTs due to large arrays. */
#if defined __MAXHISTMEM__
#undef __MAXHISTMEM__
#endif
#define __MAXHISTMEM__ 300  /* 1MB > 10^5*REAL64. 100 is sometimes too small. */
/******************************************************************************/


/******************************************************************************/
/* Will we compile the code with aggressive options (e.g., -Ofast for gfortran)?
 * Some debugging will be disabled if yes (1). Note:
 * 1. It is VALID to set __AGRESSIVE_OPTIONS__ to 0 and __DEBUGGING__ to 1 at
 * the same time.
 * 2. When compiled with aggressive options, the code may behave unexpectedly */
#if defined __AGRESSIVE_OPTIONS__
#undef __AGRESSIVE_OPTIONS__
#endif
#define __AGRESSIVE_OPTIONS__ 0
/******************************************************************************/


/******************************************************************************/
/* Do we use the STORAGE_SIZE intrinsic? (Fortran 2008) */
/* We prefer STORAGE_SIZE to C_SIZEOF, because the former is intrinsic while the
 * later requires the intrinsic module ISO_C_BINDING.
 * The flag will be removed later (in 2025?). */
#if defined __USE_STORAGE_SIZE__
#undef __USE_STORAGE_SIZE__
#endif
#define __USE_STORAGE_SIZE__ 0

#if defined __GFORTRAN__
#if __GNUC__ >= 5 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
#undef __USE_STORAGE_SIZE__
#define __USE_STORAGE_SIZE__ 1
#endif
#endif

#if defined __INTEL_COMPILER
#if __INTEL_COMPILER >= 1710
#undef __USE_STORAGE_SIZE__
#define __USE_STORAGE_SIZE__ 1
#endif
#endif

#if defined __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE >= 53
#undef __USE_STORAGE_SIZE__
#define __USE_STORAGE_SIZE__ 1
#endif
#endif

#if defined __PGI
#if __PGIC__ >= 15 && __PGIC_MINOR__ >= 4
#undef __USE_STORAGE_SIZE__
#define __USE_STORAGE_SIZE__ 1
#endif
#endif

#if defined __ibmxl__
#if __ibmxl_version__ >= 15 && __ibmxl_release__ >= 2
#undef __USE_STORAGE_SIZE__
#define __USE_STORAGE_SIZE__ 1
#endif
#endif

#if __FORTRAN_STANDARD__ < 2008
#undef __USE_STORAGE_SIZE__
#define __USE_STORAGE_SIZE__ 0
#endif
/******************************************************************************/


/******************************************************************************/
/* Do we use INT16, INT32, INT64, REAL32, REAL64, REAL128 from ISO_FORTRAN_ENV?
 * (Fortran 2008). This flag will be removed later (in 2025?). */
#if defined __USE_ISO_FORTRAN_ENV_INTREAL__
#undef __USE_ISO_FORTRAN_ENV_INTREAL__
#endif
#define __USE_ISO_FORTRAN_ENV_INTREAL__ 0

#if defined __GFORTRAN__
#if __GNUC__ >= 5 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5)
#undef __USE_ISO_FORTRAN_ENV_INTREAL__
#define __USE_ISO_FORTRAN_ENV_INTREAL__ 1
#endif
#endif

#if defined __INTEL_COMPILER
#if __INTEL_COMPILER >= 1640
#undef __USE_ISO_FORTRAN_ENV_INTREAL__
#define __USE_ISO_FORTRAN_ENV_INTREAL__ 1
#endif
#endif

#if defined __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE >= 53
#undef __USE_ISO_FORTRAN_ENV_INTREAL__
#define __USE_ISO_FORTRAN_ENV_INTREAL__ 1
#endif
#endif

#if defined __PGI
#if __PGIC__ >= 14 && __PGIC_MINOR__ >= 1
#undef __USE_ISO_FORTRAN_ENV_INTREAL__
#define __USE_ISO_FORTRAN_ENV_INTREAL__ 1
#endif
#endif

#if defined __ibmxl__
#if __ibmxl_version__ >= 14 && __ibmxl_release__ >= 1
#undef __USE_ISO_FORTRAN_ENV_INTREAL__
#define __USE_ISO_FORTRAN_ENV_INTREAL__ 1
#endif
#endif

#if __FORTRAN_STANDARD__ < 2008
#undef __USE_ISO_FORTRAN_ENV_INTREAL__
#define __USE_ISO_FORTRAN_ENV_INTREAL__ 0
#endif
/******************************************************************************/
