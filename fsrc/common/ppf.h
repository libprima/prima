#if defined __DEBUGGING__
#undef __DEBUGGING__
#endif
#define __DEBUGGING__ 1 


#if defined __FORTRAN_STANDARD__
#undef __FORTRAN_STANDARD__
#endif
#define __FORTRAN_STANDARD__ 2003 

#if defined __USE_INTRINSIC_ALGEBRA__
#undef __USE_INTRINSIC_ALGEBRA__
#endif
#define __USE_INTRINSIC_ALGEBRA__ 1
#if defined __USE_POWELL_ALGEBRA__
#undef __USE_POWELL_ALGEBRA__
#endif
#define __USE_POWELL_ALGEBRA__ 1 
#if defined __INTEGER_KIND__
#undef __INTEGER_KIND__
#endif
#define __INTEGER_KIND__ 0 
#if defined __REAL_PRECISION__
#undef __REAL_PRECISION__
#endif
#define __REAL_PRECISION__ 64 
#if defined __QP_AVAILABLE__
#undef __QP_AVAILABLE__
#endif
#define __QP_AVAILABLE__ 0  

#if __QP_AVAILABLE__ != 1 && __REAL_PRECISION__ > 64
#undef __REAL_PRECISION__
#define __REAL_PRECISION__ 64
#endif
#if defined __USE_IEEE_ARITHMETIC__
#undef __USE_IEEE_ARITHMETIC__
#endif
#define __USE_IEEE_ARITHMETIC__ 0

#if defined __GFORTRAN__
#if __REAL_PRECISION__ <= 64 && __GNUC__ >= 5
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#if defined __INTEL_COMPILER
#if __INTEL_COMPILER >= 1110
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#if defined __NAG_COMPILER_RELEASE
#if __NAG_COMPILER_RELEASE >= 50 
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#if defined __PGI
#if __PGIC__ >= 11 && __PGIC_MINOR__ >= 1
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#if defined __ibmxl__
#if __ibmxl_version__ >= 13 && __ibmxl_release__ >= 1
#undef __USE_IEEE_ARITHMETIC__
#define __USE_IEEE_ARITHMETIC__ 1
#endif
#endif

#if __FORTRAN_STANDARD__ < 2003
#undef __USE_IEEE_ARITHMETIC__ 
#define __USE_IEEE_ARITHMETIC__ 0
#endif
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
