! INFNAN is a module defining some procedures concerning INF or NAN.
!
! N.B.:
!
! 1. Here we implement the procedures for single, double, and quadruple precisions, because when we
! interface the Fortran code with other languages (e.g., MATLAB), the procedures may be invoked in
! both the Fortran code and the gateway (e.g., MEX gateway), which may use different real precisions
! (e.g., the Fortran code may use single precision, but the MEX gateway uses double by default).
!
! 2. We decide not to use IEEE_IS_NAN and IEEE_IS_FINITE provided by the intrinsic IEEE_ARITHMETIC
! available since Fortran 2003. The reason is as follows. These two procedures are supposed to
! return default logical values. However, if the code is compiled by gfortran 9.3.0 with the option
! -fdefault-integer-8, then compiler will enforce the default logical value to be 64-bit, but the
! returned kinds of IEEE_IS_NAN and IEEE_IS_FINITE will not be changed accordingly, and they will
! remain 32-bit if that is the default logical kind without -fdefault-integer-8. Consequently, the
! returned kinds of IEEE_IS_NAN and IEEE_IS_FINITE may actually differ from the default logical kind
! due to this compiler option and hence violate the Fortran standard! This is fatal, because a
! piece of perfectly standard-compliant may fail to be compiled due to type mismatches. It is
! similar with ifort 2021.2.0 and nagfor 7.0. See the following address for more discussions:
! https://stackoverflow.com/questions/69060408/bug-matlab-mex-changes-the-kind-of-the-default-logical
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Sunday, September 05, 2021 PM10:30:56


#include "ppf.h"

module infnan_mod

implicit none
private
public :: is_nan, is_finite, is_posinf, is_neginf, is_inf

#if __QP_AVAILABLE__ == 1

interface is_nan
    module procedure is_nan_sp, is_nan_dp, is_nan_qp
end interface is_nan

interface is_finite
    module procedure is_finite_sp, is_finite_dp, is_finite_qp
end interface is_finite

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

interface is_nan
    module procedure is_nan_sp, is_nan_dp
end interface is_nan

interface is_finite
    module procedure is_finite_sp, is_finite_dp
end interface is_finite

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
! nan, or finite.


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
