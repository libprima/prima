! WARNERROR is a module defining some procedures concerning warnings and
! errors.

#include "ppf.h"

      module warnerror_mod

      implicit none
      private
      public :: errstop, backtr

      contains
          
      subroutine errstop(srname, message)
      character(len=*), intent(in) :: srname
      character(len=*), intent(in) :: message 

      call backtr()

      print *
      print '(1x, 5A)', 'Error: ',  srname(1 : len_trim(srname)),       &
     & ': ', message(1 : len_trim(message)), '.'
      print *

      stop

      end subroutine errstop

      subroutine backtr()
      ! BACKTR calls a compiler-dependent intrinsic to show a backtrace
      ! if we are in the debuge mode, i.e., __DEBUG_MODE__ == 1.
      ! N.B.: The intrisic depends on compiler and does not exist in all
      ! compilers. Indeed, it is not standard-conforming. Therefore,
      ! compilers may warn that a non-standard intrinsic is in use.
      implicit none
#if __DEBUG_MODE__ == 1 

#ifdef __GFORTRAN__
      call backtrace()
#endif

#ifdef __INTEL_COMPILER
      call tracebackqq()
#endif

#endif
      end subroutine backtr

      end module warnerror_mod
