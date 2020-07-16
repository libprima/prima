/*************************************************************************/
/* 
ppf.h defines the following preprocessing flags (the first value is default). 

__DEBUG_MODE__              debug or not: 0, 1
__INTEGER_KIND__            the integer kind to be used: 0, 32, 16, 64
__REAL_PRECISION__          the real precision to be used: 64, 32, 128, 0 
__FORTRAN_STANDARD__        Fortran standard to follow: 95, 2003, 2008
__IMPROVE_POWELL_CODE__     improve Powell's code or not: 1, 0
__USE_IEEE_ARITHMETIC__     use the IEEE_ARITHMETIC intrinsic or not: 1, 0
__USE_INTRINSIC_ALGEBRA__   use intrinsic procedures like matmul or not: 1, 0

You may change these flags, but make sure that your compiler is supportive
when changing __INTEGER_KIND__, __REAL_PRECISION__, __FORTRAN_STANDARD__,
and __USE_INTRINSIC_ALGEBRA__.

Why not define these flags as parameters in the Fortran code, e.g.,

logical, parameter :: __DEBUG_MODE__ = .false. ?

Such a definition will work for __DEBUG_MODE__, but not for the flags that
depend on the compiler, for instance, __FORTRAN_STANDARD__.
*/
/*************************************************************************/


/*************************************************************************/
/* Are we debugging? */
#ifdef __DEBUG_MODE__
#undef __DEBUG_MODE__
#endif
#define __DEBUG_MODE__ 1 
/*************************************************************************/


/*************************************************************************/
/* Which integer kind to use? 
 * 0 = default INtEGER, 16 = INTEGER*2, 32 = INTEGER*4, 64 = INTEGER*8. 
 * Make sure that your compiler supports the selected kind. */
#ifdef __INTEGER_KIND__
#undef __INTEGER_KIND__
#endif
#define __INTEGER_KIND__ 16 
/*************************************************************************/


/*************************************************************************/
/* Which real kind to use? 
 * 0 = default REAL, 32 = REAL*4, 64 = REAL*8, 128 = REAL*16.
 * Make sure that your compiler supports the selected kind. 
 * Note: The default REAL (i.e., 0) is the single-precision REAL. */
#ifdef __REAL_PRECISION__
#undef __REAL_PRECISION__
#endif
#define __REAL_PRECISION__ 64 

/* Is quad precision available on this platform (compiler, hardware ...)? */
/* Note:
 * 1. Not all platforms support REAL128. For example, pgfortran 19 does not.
 * 2. It is not guaranteed that REAL128 has a wider range than REAL64. For
 *    example, REAL128 of nagfor 7.0 has a range of 291, while REAL64
 *    has a range of 307. 
 * 3. It is rarely a good idea to use REAL128 as the working precision,
 *    which is probably inefficient and unnecessary. */
#ifdef __QP_AVAILABLE__
#undef __QP_AVAILABLE__
#endif
#define __QP_AVAILABLE__ 0

/* Revise __REAL_PRECISION__ according to __QP_AVAILABLE__ . */
#if __QP_AVAILABLE__ != 1 && __REAL_PRECISION__ > 64
#undef __REAL_PRECISION__
#define __REAL_PRECISION__ 64
#endif
/*************************************************************************/


/*************************************************************************/
/* Which Fortran standard do we follow? */
/* We aim to be compatible with Fortran 95, 2003 and 2008. 
 * Make sure that your compiler supports the selected standard. */
#ifdef __FORTRAN_STANDARD__
#undef __FORTRAN_STANDARD__
#endif
/*#define __FORTRAN_STANDARD__ 95 */
#define __FORTRAN_STANDARD__ 2003
/*#define __FORTRAN_STANDARD__ 2008 */

/* Revise __FORTRAN_STANDARD__ according to the version of the compiler. */
/* Of course, we cannot exhaust all the compilers. */
/*******************************************************/
#if __FORTRAN_STANDARD__ > 2003

#ifdef __GFORTRAN__
#if __GNUC__ < 5  /* gfortran 5.0 supports most features of F2003/08 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95
#endif
#endif

#ifdef __INTEL_COMPILER
#if __INTEL_COMPILER < 1800 /* ifort 18.0 fully supports F2008 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 2003 
#endif
#endif

#ifdef __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE < 62  /* nagfor 6.2 supports most of F2008 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 2003 
#endif
#endif

#ifdef __PGI
/* pgifortran 11 fully supports F2003; support for F2008 is increasing */
#if __PGIC__ < 11  
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#ifdef __G95__
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95  /* g95 only supports F2003/08 partially */
#endif

#endif
/*******************************************************/

/*******************************************************/
#if __FORTRAN_STANDARD__ > 95

#ifdef __GFORTRAN__
#if __GNUC__ < 5  /* gfortran 5.0 supports most of F2003/08 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#ifdef __INTEL_COMPILER
#if __INTEL_COMPILER < 1600  /* ifort 16.0 fully supports F2003 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#ifdef __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE < 61  /* nagfor 6.1 fully supports F2003 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#ifdef __PGI
#if __PGIC__ < 11  /* pgifortran 11 fully supports F2003 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#ifdef __G95__
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 /* g95 only supports F2003/08 partially */
#endif

#endif
/*******************************************************/
/*************************************************************************/


/*************************************************************************/
/* Do we improve Powell's code? */
/* The improvements do not change Powell's algorithms but modify the 
 * implementation of some algebraic calculations. The improved code may
 * not produce exactly the same results as Powell's code due to properties
 * of floating-point arithmetic, e.g., the non-associativity of
 * floating-point addition and multiplication. The improvements remove
 * some bugs. */
#ifdef __IMPROVE_POWELL_CODE__
#undef __IMPROVE_POWELL_CODE__
#endif
#define __IMPROVE_POWELL_CODE__ 0
/*************************************************************************/


/*************************************************************************/
/* Do we use IEEE_ARITHMETIC? */
/* Make sure that your compiler supports IEEE_ARITHMETIC if you set this
 * value to 1. */
#ifdef __USE_IEEE_ARITHMETIC__      
#undef __USE_IEEE_ARITHMETIC__
#endif
#if __FORTRAN_STANDARD__ >= 2003 
/* IEEE_ARITHMETIC is available starting from Fortran 2003. */
#define __USE_IEEE_ARITHMETIC__ 1 
#else
#define __USE_IEEE_ARITHMETIC__ 0 
#endif
/* As of gfortran 5.5, it seems that the IEEE_ARITHMETIC of gfortran does 
 * not support REAL128. */
#if __REAL_PRECISION__ > 64
#ifdef __GNUC__
#undef __USE_IEEE_ARITHMETIC__ 
#define __USE_IEEE_ARITHMETIC__ 0
#endif
#endif
/*************************************************************************/


/*************************************************************************/
/* Do we use the intrinsic algebra procedures (e.g., matmul)? */
/* If no, we use the procedures implemented in lina.F. */
#ifdef __USE_INTRINSIC_ALGEBRA__
#undef __USE_INTRINSIC_ALGEBRA__
#endif
#define __USE_INTRINSIC_ALGEBRA__ 0 
/* We do not use intrinsic algebra procedures in debug mode. Instead, we
 * use our own implementation of these procedures in lina.F. */
#if __DEBUG_MODE__ == 1
#undef __USE_INTRINSIC_ALGEBRA__
#define __USE_INTRINSIC_ALGEBRA__ 0
#endif
/*************************************************************************/
