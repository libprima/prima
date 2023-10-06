#include "ppf.h"

module huge_mod
!--------------------------------------------------------------------------------------------------!
! This module provides a function that returns HUGE(X). See infnan.f90 for more comments.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Friday, October 06, 2023 PM07:30:09
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: huge_value

#if PRIMA_QP_AVAILABLE == 1

interface huge_value
    module procedure huge_value_sp, huge_value_dp, huge_value_qp
end interface huge_value

#else

interface huge_value
    module procedure huge_value_sp, huge_value_dp
end interface huge_value

#endif


contains


pure elemental function huge_value_sp(x) result(y)
use, non_intrinsic :: consts_mod, only : SP
implicit none
real(SP), intent(in) :: x
real(SP) :: y
y = huge(x)
end function huge_value_sp

pure elemental function huge_value_dp(x) result(y)
use, non_intrinsic :: consts_mod, only : DP
implicit none
real(DP), intent(in) :: x
real(DP) :: y
y = huge(x)
end function huge_value_dp

#if PRIMA_QP_AVAILABLE == 1
pure elemental function huge_value_qp(x) result(y)
use, non_intrinsic :: consts_mod, only : QP
implicit none
real(QP), intent(in) :: x
real(QP) :: y
y = huge(x)
end function huge_value_qp
#endif

end module huge_mod
