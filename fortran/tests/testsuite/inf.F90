#include "ppf.h"

module inf_mod
!--------------------------------------------------------------------------------------------------!
! This is the replacement of fortran/common/inf.F90 for testing what if IS_FINITE always returns
! TRUE, while is_inf, is_posinf, is_neginf, and is_nan always return FALSE.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Friday, October 06, 2023 PM08:29:32
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: huge_mod, only : huge_value
implicit none
private
public :: is_finite, is_posinf, is_neginf, is_inf

#if PRIMA_QP_AVAILABLE == 1

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


pure elemental function is_finite_sp(x) result(y)
use, non_intrinsic :: consts_mod, only : SP
implicit none
real(SP), intent(in) :: x
logical :: y
y = .true.
end function is_finite_sp

pure elemental function is_finite_dp(x) result(y)
use, non_intrinsic :: consts_mod, only : DP
implicit none
real(DP), intent(in) :: x
logical :: y
y = .true.
end function is_finite_dp

pure elemental function is_posinf_sp(x) result(y)
use, non_intrinsic :: consts_mod, only : SP
implicit none
real(SP), intent(in) :: x
logical :: y
y = .false.
end function is_posinf_sp

pure elemental function is_posinf_dp(x) result(y)
use, non_intrinsic :: consts_mod, only : DP
implicit none
real(DP), intent(in) :: x
logical :: y
y = .false.
end function is_posinf_dp

pure elemental function is_neginf_sp(x) result(y)
use, non_intrinsic :: consts_mod, only : SP
implicit none
real(SP), intent(in) :: x
logical :: y
y = .false.
end function is_neginf_sp

pure elemental function is_neginf_dp(x) result(y)
use, non_intrinsic :: consts_mod, only : DP
implicit none
real(DP), intent(in) :: x
logical :: y
y = .false.
end function is_neginf_dp

pure elemental function is_inf_sp(x) result(y)
use, non_intrinsic :: consts_mod, only : SP
implicit none
real(SP), intent(in) :: x
logical :: y
y = .false.
end function is_inf_sp

pure elemental function is_inf_dp(x) result(y)
use, non_intrinsic :: consts_mod, only : DP
implicit none
real(DP), intent(in) :: x
logical :: y
y = .false.
end function is_inf_dp


#if PRIMA_QP_AVAILABLE == 1

pure elemental function is_finite_qp(x) result(y)
use, non_intrinsic :: consts_mod, only : QP
implicit none
real(QP), intent(in) :: x
logical :: y
y = .true.
end function is_finite_qp

pure elemental function is_posinf_qp(x) result(y)
use, non_intrinsic :: consts_mod, only : QP
implicit none
real(QP), intent(in) :: x
logical :: y
y = .false.
end function is_posinf_qp

pure elemental function is_neginf_qp(x) result(y)
use, non_intrinsic :: consts_mod, only : QP
implicit none
real(QP), intent(in) :: x
logical :: y
y = .false.
end function is_neginf_qp

pure elemental function is_inf_qp(x) result(y)
use, non_intrinsic :: consts_mod, only : QP
implicit none
real(QP), intent(in) :: x
logical :: y
y = .false.
end function is_inf_qp

#endif


end module inf_mod
