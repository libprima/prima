#include "ppf.h"
program testsize
     !use, intrinsic :: iso_c_binding, only : c_sizeof
#if __FORTRAN_STANDARD__ >= 2008
      use, intrinsic :: iso_fortran_env, only : INT16, INT32, INT64, REAL32, REAL64, REAL128 
      implicit none
#else
      implicit none
    integer, parameter :: REAL32 = kind(0.0) 
    integer, parameter :: REAL64 = kind(0.0D0)
    integer, parameter :: REAL128 = selected_real_kind(p = 30)
#endif
      real(REAL32) :: xs
      real(REAL64) :: xd
      real(REAL128) :: xq
      print *, storage_size(xs)!, c_sizeof(xs)
      print *, storage_size(xd)!, c_sizeof(xd)
      print *, storage_size(xq)!, c_sizeof(xq)
end program testsize
