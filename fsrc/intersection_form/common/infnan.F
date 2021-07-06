!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of infnan.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 06-Jul-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! INFNAN is a module defining some procedures concerning INF or NAN.
!
! Here we implement the procedures for single, double, and quadruple
! precisions, because when we interface the Fortran code with other
! languages (e.g., MATLAB), the procedures may be invoked in both the
! Fortran code and the gateway (e.g., MEX gateway), which may use
! different real precisions (e.g., the Fortran code may use single
! precision, but MEX gateway uses double precision by default).
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Sunday, May 23, 2021 AM10:59:55


#include "ppf.h"

      module infnan_mod

#if __USE_IEEE_ARITHMETIC__ != 0
      use, intrinsic :: ieee_arithmetic, only : is_nan => ieee_is_nan, i&
     &s_finite => ieee_is_finite
#endif

      implicit none
      private
      public :: is_nan, is_finite, is_posinf, is_neginf, is_inf

#if __QP_AVAILABLE__ == 1

#if __USE_IEEE_ARITHMETIC__ == 0
      interface is_nan
          module procedure is_nan_sp, is_nan_dp, is_nan_qp
      end interface is_nan
      interface is_finite
          module procedure is_finite_sp, is_finite_dp, is_finite_qp
      end interface is_finite
#endif

      interface is_posinf
          module procedure is_posinf_sp, is_posinf_dp, is_posinf_qp
      end interface is_posinf
      interface is_neginf
          module procedure is_neginf_sp, is_neginf_dp, is_neginf_qp
      end interface is_neginf
      interface is_inf
          module procedure is_inf_sp, is_inf_dp, is_inf_qp
      end interface is_inf

#else

#if __USE_IEEE_ARITHMETIC__ == 0
      interface is_nan
          module procedure is_nan_sp, is_nan_dp
      end interface is_nan
      interface is_finite
          module procedure is_finite_sp, is_finite_dp
      end interface is_finite
#endif

      interface is_posinf
          module procedure is_posinf_sp, is_posinf_dp
      end interface is_posinf
      interface is_neginf
          module procedure is_neginf_sp, is_neginf_dp
      end interface is_neginf
      interface is_inf
          module procedure is_inf_sp, is_inf_dp
      end interface is_inf
#endif


      contains
! The following functions check whether a real number x is infinite,
! nan, or finite. Starting from Fortran 2003, the intrincic
! ieee_arithmetic module provides ieee_is_nan() and ieee_is_finite().


#if __USE_IEEE_ARITHMETIC__ == 0
      pure elemental function is_nan_sp(x) result(y)
      use consts_mod, only : SP
      implicit none
      real(SP), intent(in) :: x
      logical :: y
      y = (.not. (x >= x))
! This implementation may be more costly than y = (x /= x).
! However, if we defined is_nan = (x /= x), gfortran will complain
! about inequality comparison between floating-point numbers.
      end function is_nan_sp

      pure elemental function is_nan_dp(x) result(y)
      use consts_mod, only : DP
      implicit none
      real(DP), intent(in) :: x
      logical :: y
      y = (.not. (x >= x))
      end function is_nan_dp

      pure elemental function is_finite_sp(x) result(y)
      use consts_mod, only : SP
      implicit none
      real(SP), intent(in) :: x
      logical :: y
      y = (abs(x) <= huge(x))
      end function is_finite_sp

      pure elemental function is_finite_dp(x) result(y)
      use consts_mod, only : DP
      implicit none
      real(DP), intent(in) :: x
      logical :: y
      y = (abs(x) <= huge(x))
      end function is_finite_dp

#endif

      pure elemental function is_posinf_sp(x) result(y)
      use consts_mod, only : SP
      implicit none
      real(SP), intent(in) :: x
      logical :: y
      y = (x > huge(x))
      end function is_posinf_sp

      pure elemental function is_posinf_dp(x) result(y)
      use consts_mod, only : DP
      implicit none
      real(DP), intent(in) :: x
      logical :: y
      y = (x > huge(x))
      end function is_posinf_dp

      pure elemental function is_neginf_sp(x) result(y)
      use consts_mod, only : SP
      implicit none
      real(SP), intent(in) :: x
      logical :: y
      y = (-x > huge(x))
      end function is_neginf_sp

      pure elemental function is_neginf_dp(x) result(y)
      use consts_mod, only : DP
      implicit none
      real(DP), intent(in) :: x
      logical :: y
      y = (-x > huge(x))
      end function is_neginf_dp

      pure elemental function is_inf_sp(x) result(y)
      use consts_mod, only : SP
      implicit none
      real(SP), intent(in) :: x
      logical :: y
      y = (abs(x) > huge(x))
      end function is_inf_sp

      pure elemental function is_inf_dp(x) result(y)
      use consts_mod, only : DP
      implicit none
      real(DP), intent(in) :: x
      logical :: y
      y = (abs(x) > huge(x))
      end function is_inf_dp

#if __QP_AVAILABLE__ == 1
#if __USE_IEEE_ARITHMETIC__ == 0
      pure elemental function is_nan_qp(x) result(y)
      use consts_mod, only : QP
      implicit none
      real(QP), intent(in) :: x
      logical :: y
      y = (.not. (x >= x))
      end function is_nan_qp

      pure elemental function is_finite_qp(x) result(y)
      use consts_mod, only : QP
      implicit none
      real(QP), intent(in) :: x
      logical :: y
      y = (abs(x) <= huge(x))
      end function is_finite_qp
#endif

      pure elemental function is_neginf_qp(x) result(y)
      use consts_mod, only : QP
      implicit none
      real(QP), intent(in) :: x
      logical :: y
      y = (-x > huge(x))
      end function is_neginf_qp

      pure elemental function is_posinf_qp(x) result(y)
      use consts_mod, only : QP
      implicit none
      real(QP), intent(in) :: x
      logical :: y
      y = (x > huge(x))
      end function is_posinf_qp

      pure elemental function is_inf_qp(x) result(y)
      use consts_mod, only : QP
      implicit none
      real(QP), intent(in) :: x
      logical :: y
      y = (abs(x) > huge(x))
      end function is_inf_qp
#endif


      end module infnan_mod