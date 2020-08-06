/*************************************************************************/
/* 
ppf.h defines the following preprocessing flags (the first value is default). 

__DEBUGGING__               debug or not: 0, 1
__USE_INTRINSIC_ALGEBRA__   use intrinsic procedures like matmul or not: 0, 1 
__USE_POWELL_ALGEBRA__      use Powell's linear algebra procedures or not: 0, 1 
__INTEGER_KIND__            the integer kind to be used: 0, 32, 16, 64
__REAL_PRECISION__          the real precision to be used: 64, 32, 128, 0 
__FORTRAN_STANDARD__        Fortran standard to follow: 95, 2003, 2008
__USE_IEEE_ARITHMETIC__     use the IEEE_ARITHMETIC intrinsic or not: 1, 0

You may change these flags, but make sure that your compiler is supportive
when changing __INTEGER_KIND__, __REAL_PRECISION__, __FORTRAN_STANDARD__,
and __USE_INTRINSIC_ALGEBRA__.

Why not define these flags as parameters in the Fortran code, e.g.,

logical, parameter :: __DEBUGGING__ = .false. ?

Such a definition will work for __DEBUGGING__, but not for the flags that
depend on the compiler, for instance, __FORTRAN_STANDARD__.
*/
/*************************************************************************/


/*************************************************************************/
/* Are we debugging? */
#if defined __DEBUGGING__
#undef __DEBUGGING__
#endif
#define __DEBUGGING__ 1 
/*************************************************************************/


/*************************************************************************/
/* Do we use the intrinsic algebra procedures (e.g., matmul)? */
/* If no, we use the procedures implemented in lina.F. */
#if defined __USE_INTRINSIC_ALGEBRA__
#undef __USE_INTRINSIC_ALGEBRA__
#endif
#define __USE_INTRINSIC_ALGEBRA__ 0 
/*************************************************************************/


/*************************************************************************/
/* Do we use Powell's linear algebra procedures? */
/* If not, the implementation of some algebraic calculations will be
 * modified, mainly by replacing loops with matrix-vector operations.
 * This does not change Powell's algorithms, but it may not produce
 * exactly the same results as Powell's code due to properties of 
 * floating-point arithmetic, e.g., the non-associativity of floating-point
 * addition and multiplication. */
#if defined __USE_POWELL_ALGEBRA__
#undef __USE_POWELL_ALGEBRA__
#endif
#define __USE_POWELL_ALGEBRA__ 1 
/*************************************************************************/


/*************************************************************************/
/* Which integer kind to use? 
 * 0 = default INTEGER, 16 = INTEGER*2, 32 = INTEGER*4, 64 = INTEGER*8. 
 * Make sure that your compiler supports the selected kind. */
#if defined __INTEGER_KIND__
#undef __INTEGER_KIND__
#endif
#define __INTEGER_KIND__ 16 
/*************************************************************************/


/*************************************************************************/
/* Which real kind to use? 
 * 0 = default REAL, 32 = REAL*4, 64 = REAL*8, 128 = REAL*16.
 * Make sure that your compiler supports the selected kind. 
 * Note: The default REAL (i.e., 0) is the single-precision REAL. */
#if defined __REAL_PRECISION__
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
#if defined __QP_AVAILABLE__
#undef __QP_AVAILABLE__
#endif
/* Change the following line to set __QP_AVAILABLE__ to 1 if REAL128 is
 * available and if you intend to use it. 
 * Note:
 * 1. Do NOT change is unless you are really sure.
 * 2. If the code is terfaced with MATLAB, then you still need to modify
 *    mexapi.F to set __USE_QP__ to 1. */
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
#if defined __FORTRAN_STANDARD__
#undef __FORTRAN_STANDARD__
#endif
/*#define __FORTRAN_STANDARD__ 95 */
#define __FORTRAN_STANDARD__ 2003
/*#define __FORTRAN_STANDARD__ 2008 */

/* Revise __FORTRAN_STANDARD__ according to the version of the compiler. */
/* Of course, we cannot exhaust all the compilers. */
/*******************************************************/
#if __FORTRAN_STANDARD__ > 2003

#if defined __GFORTRAN__
#if __GNUC__ < 5  /* gfortran 5.0 supports most features of F2003/08 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95
#endif
#endif

#if defined __INTEL_COMPILER
#if __INTEL_COMPILER < 1800  /* ifort 18.0 fully supports F2008 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 2003 
#endif
#endif

#if defined __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE < 62  /* nagfor 6.2 supports most of F2008 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 2003 
#endif
#endif

#if defined __PGI
/* pgifortran 11 fully supports F2003; support for F2008 is increasing */
#if __PGIC__ < 11  
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#if defined __G95__
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95  /* g95 only supports F2003/08 partially */
#endif

#endif
/*******************************************************/

/*******************************************************/
#if __FORTRAN_STANDARD__ > 95

#if defined __GFORTRAN__
#if __GNUC__ < 5  /* gfortran 5.0 supports most of F2003/08 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#if defined __INTEL_COMPILER
#if __INTEL_COMPILER < 1600  /* ifort 16.0 fully supports F2003 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#if defined __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE < 61  /* nagfor 6.1 fully supports F2003 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#if defined __PGI
#if __PGIC__ < 11  /* pgifortran 11 fully supports F2003 */
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 
#endif
#endif

#if defined __G95__
#undef __FORTRAN_STANDARD__
#define __FORTRAN_STANDARD__ 95 /* g95 only supports F2003/08 partially */
#endif

#endif
/*******************************************************/
/*************************************************************************/


/*************************************************************************/
/* Do we use IEEE_ARITHMETIC? */
/* Make sure that your compiler supports IEEE_ARITHMETIC if you set this
 * value to 1. */
#if defined __USE_IEEE_ARITHMETIC__      
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
#if defined __GNUC__
#undef __USE_IEEE_ARITHMETIC__ 
#define __USE_IEEE_ARITHMETIC__ 0
#endif
#endif
/*************************************************************************/
