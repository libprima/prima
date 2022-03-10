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
 * __FORTRAN_STANDARD__        which Fortran standard to follow: 2003, 2008, 2018
 * __USE_POWELL_ALGEBRA__      use Powell's linear algebra procedures or not: 1, 0
 * __USE_INTRINSIC_ALGEBRA__   use intrinsic procedures like matmul or not: 0, 1
 * __INTEGER_KIND__            the integer kind to be used: 0, 32, 64, 16
 * __REAL_PRECISION__          the real precision to be used: 64, 32, 128, 0
 * __USE_STORAGE_SIZE__        use the STORAGE_SIZE intrinsic or not: 0, 1
 * __USE_ISO_FORTRAN_ENV_INTREAL__ use INT32 etc in ISO_FORTRAN_ENV or not: 0, 1
 *
 * N.B.:
 *
 * 0. USE THE DEFAULT IF UNSURE.
 *
 * 1. When __USE_POWELL_ALGEBRA__ == 1, the released version will produce EXACTLY
 * the same results as Powell's Fortran 77 code (necessarily, this means that the
 * code performs EXACTLY the same calculations as Powell's code, or else rounding
 * errors will lead to different computed results, the difference being sometimes
 * significant for nonconvex problems).
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
#if defined RELEASED__
#undef RELEASED__
#endif
#define RELEASED__ 0
/******************************************************************************/


/******************************************************************************/
/* Are we debugging?
 * __RELEASED__ == 1 and __DEBUGGING__ == 1 do not conflict. Users may debug. */
#if defined __DEBUGGING__
#undef __DEBUGGING__
#endif
#define __DEBUGGING__ 1
/******************************************************************************/


/******************************************************************************/
/* Which Fortran standard to follow? */
#if defined __FORTRAN_STANDARD__
#undef __FORTRAN_STANDARD__
#endif
#define __FORTRAN_STANDARD__ 2003 /* Will be default to 2008 later (in 2025?).*/
/******************************************************************************/


/******************************************************************************/
/* Do we use Powell's linear algebra procedures? */
/* If not, some basic algebraic procedures will be implemented with matrix/vector
 * operations instead of loops. This does not change Powell's algorithms, but
 * it may not produce exactly the same results as Powell's code due to properties
 * of floating-point arithmetic, e.g., the non-associativity of floating-point
 * addition and multiplication. */
#if defined __USE_POWELL_ALGEBRA__
#undef __USE_POWELL_ALGEBRA__
#endif
#define __USE_POWELL_ALGEBRA__ 1
/******************************************************************************/


/******************************************************************************/
/* Do we use the intrinsic algebra procedures (e.g., matmul, dot_product)? */
/* If no, we use the procedures implemented in linalg.F. */
/* When __USE_INTRINSIC_ALGEBRA__ == 1, the code may not produce exactly the
 * same results as Powell's code, because the intrinsic matmul behaves
 * differently from a naive triple loop due to finite-precision arithmetic.
 * The difference has been observed on matprod22 and matprod12. The second case
 * occurred on Oct. 11, 2021 in the trust-region subproblem solver of COBYLA, and
 * it took enormous time to find out that Powell's code and the modernized code
 * behaved differently due to matmul and matprod12 when calculating RESMAX (in
 * Powell's code) and CSTRV (in the modernized code) when stage 2 starts. */
#if defined __USE_INTRINSIC_ALGEBRA__
#undef __USE_INTRINSIC_ALGEBRA__
#endif
#define __USE_INTRINSIC_ALGEBRA__ 0
/******************************************************************************/


/******************************************************************************/
/* __USE_POWELL_ALGEBRA__ == 1 and __USE_INTRINSIC_ALGEBRA__ == 1 do NOT conflict.
 * However, to make sure that the code produces exactly the same results as
 * Powell's code when __USE_POWELL_ALGEBRA__ == 1, we impose the following.*/
#if __USE_POWELL_ALGEBRA__ == 1 && __RELEASED__ == 1
#undef __USE_INTRINSIC_ALGEBRA__
#define __USE_INTRINSIC_ALGEBRA__ 0
#endif
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
