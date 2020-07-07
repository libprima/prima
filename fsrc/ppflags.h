/* PPFLAGS defines some flags for preprocessing. */

/* Are we debugging? */
#ifdef __DEBUG_MODE__
#undef __DEBUG_MODE__
#endif
#define __DEBUG_MODE__ 1


/* Which Fortran standard do we follow? */
/* We aim to be compatible with Fortran 95, 2003 and 2008. However, the
 * compiler must support the selected standard. Otherwiswe, it will not
 * work. */
#ifdef __FORTRAN_STANDARD__
#undef __FORTRAN_STANDARD__
#endif
#define __FORTRAN_STANDARD__ 2008 
/*#define __FORTRAN_STANDARD__ 2003 */
/*#define __FORTRAN_STANDARD__ 2008 */ 


/* Revise __FORTRAN_STANDARD__ according to the version of the compiler. */
/* Of course, we cannot exhaust all the compilers. */
/*************************************************************************/
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
/* pgifortran 11 fully supports F2003; support for F2008 is increasing. */
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

/*************************************************************************/
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
/*************************************************************************/


/* Do we intend to interface the code with MATLAB via MEX? */
/* If yes, then the integer type used in the Fortran code must be the
 * default integer. Otherwise, it may not be compatible with MATLAB MEX
 * and lead to Segmentation Faults. */
#ifdef __INTEFACE_WITH_MATLAB__
#undef __INTEFACE_WITH_MATLAB__
#endif
#define __INTEFACE_WITH_MATLAB__ 1


/* Do we improve Powell's code? */
/* The improvements do not change Powell's algorithms but modify some
 * implementations of algebraic calculations. The improved code may not
 * produce exactly the same results as Powell's code due to properties
 * of floating-point arithmetic (e.g., the non-associativity of
 * floating-point addition and multiplication. Some bugs will also be
 * removed by the improvements. */
#ifdef __IMPROVE_POWELL_CODE__
#undef __IMPROVE_POWELL_CODE__
#endif
#define __IMPROVE_POWELL_CODE__ 0


/* Do we use the intrinsic algebra procedures (e.g., matmul)? */
/* If no, we use the procedures implemented in lina.F. */
#ifdef __USE_INTRINSIC_ALGEBRA__
#undef __USE_INTRINSIC_ALGEBRA__
#endif
#define __USE_INTRINSIC_ALGEBRA__ 0 
/* We do not use intrinsic algebra procedures in debug mode. Instead, we
 * use our own implementation of these procedures. */
#if __DEBUG_MODE__ == 1
#undef __USE_INTRINSIC_ALGEBRA__
#define __USE_INTRINSIC_ALGEBRA__ 0
#endif


/* Do we use IEEE_ARITHMETIC? */
#ifdef __USE_IEEE_ARITHMETIC__      
#undef __USE_IEEE_ARITHMETIC__
#endif
#define __USE_IEEE_ARITHMETIC__ 0
/*************************************************************************/
/* Decide whether IEEE_ARITHMETIC is supported the compiler. */
/* Of course, we cannot exhaust all the compilers. */
#if __FORTRAN_STANDARD__ >= 2003 
/* IEEE_ARITHMETIC is available starting from Fortran 2003. */

#ifdef __GFORTRAN__
#if __GNUC__ > 4
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#ifdef __INTEL_COMPILER
#if __INTEL_COMPILER > 1100
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#ifdef __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE > 40
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#ifdef __G95__
#if __G95_MINOR__ > 90
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#endif
/*************************************************************************/
