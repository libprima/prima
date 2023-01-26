#include "ppf.h"

module infnan_mod
!--------------------------------------------------------------------------------------------------!
! This module provides functions that check whether a real number X is infinite, NaN, or finite.
!
! N.B.:
!
! 1. We implement all the procedures for single, double, and quadruple precisions (when available).
! When we interface the Fortran code with other languages (e.g., MATLAB), the procedures may be
! invoked in both the Fortran code and the gateway (e.g., MEX gateway), which may use different real
! precisions (e.g., the Fortran code may use single, but the MEX gateway uses double by default).
!
! 2. We decide not to use IEEE_IS_NAN and IEEE_IS_FINITE provided by the intrinsic IEEE_ARITHMETIC
! available since Fortran 2003. The reason is as follows. The Fortran standards require these two
! procedures to return default logical values. However, if the code is compiled by gfortran 9.3.0
! with the option -fdefault-integer-8 (which is adopted by MEX and cannot be changed easily), then
! the compiler will enforce the default logical value to be 64-bit, but the returned kinds of
! IEEE_IS_NAN and IEEE_IS_FINITE will not be changed accordingly, and they will remain 32-bit if
! that is the default logical kind. Therefore, the returned kinds of IEEE_IS_NAN and IEEE_IS_FINITE
! may actually differ from the default logical kind due to this compiler option and hence violate
! the Fortran standard! This is fatal, because a piece of perfectly standard-compliant code may fail
! to be compiled due to type mismatches. It is similar with ifort 2021.2.0 and nagfor 7.0. See more
! discussions at https://stackoverflow.com/questions/69060408.
!
! 3. The functions aim to work even when compilers are invoked with aggressive optimization flags,
! such as `gfortran -Ofast`.
!
! 4. There are many ways to implement functions like IS_NAN. However, not all of them work with
! aggressive optimization flags. For example, for gfortran 9.3.0, the IEEE_IS_NAN included in
! IEEE_ARITHMETIC does not work with `gfortran -Ofast`. Another example, when X is NaN, (X == X) and
! (X >= X) are evaluated as TRUE by Flang 7.1.0 and nvfortran 21.3-0, even if they are invoked
! without any explicit optimization flag. See the following for discussions
! https://stackoverflow.com/questions/15944614
!
! 5. The most naive implementation for IS_NAN is (X /= X). However, compilers (e.g., gfortran) may
! complain about inequality comparison between floating-point numbers. In addition, it is likely to
! fail when compliers are invoked with aggressive optimization flags.
!
! 6. The implementation below is totally empirical, in the sense that I have not studied in-depth
! what the aggressive optimization flags really do, but only made some tests and found the
! implementation that worked correctly. The story may change when compilers are changed/updated.
!
! 7. N.B.: Do NOT change the functions without thorough testing. Their implementations are delicate.
! For example, when compilers are invoked with aggressive optimization flags,
! (X <= HUGE(X) .AND. X >= -HUGE(X)) may differ from (ABS(X) <= HUGE(X)) ,
! (X > HUGE(X) .OR. X < -HUGE(X)) may differ from (ABS(X) > HUGE(X)) , and
! (ABS(X) > HUGE(X) .AND. X > 0) may differ from (X > HUGE(X)) .
!
! 8. IS_NAN must be implemented in a file separated from is_inf and is_finite (a separated module is
! not enough). Otherwise, IS_NAN may not work with some compilers invoked with aggressive
! optimization flags e.g., ifx -fast with ifx 2022.1.0 or flang -Ofast with flang 15.0.3.
!
! 9. Even though the functions involve invocation of ABS and HUGE, their performance (in terms of
! CPU time) turns out comparable to or even better than the functions in IEEE_ARITHMETIC.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Tuesday, January 24, 2023 PM04:17:14
!--------------------------------------------------------------------------------------------------!

use inf_mod, only : is_finite, is_inf, is_posinf, is_neginf
implicit none
private
public :: is_finite, is_posinf, is_neginf, is_inf, is_nan

#if __QP_AVAILABLE__ == 1

interface is_nan
    module procedure is_nan_sp, is_nan_dp, is_nan_qp
end interface is_nan

#else

interface is_nan
    module procedure is_nan_sp, is_nan_dp
end interface is_nan

#endif


contains


pure elemental function is_nan_sp(x) result(y)
use, non_intrinsic :: consts_mod, only : SP
implicit none
real(SP), intent(in) :: x
logical :: y
!y = (.not. (x <= huge(x) .and. x >= -huge(x))) .and. (.not. abs(x) > huge(x))  ! Does not always work
y = (.not. is_finite(x)) .and. (.not. is_inf(x))
end function is_nan_sp

pure elemental function is_nan_dp(x) result(y)
use, non_intrinsic :: consts_mod, only : DP
implicit none
real(DP), intent(in) :: x
logical :: y
!y = (.not. (x <= huge(x) .and. x >= -huge(x))) .and. (.not. abs(x) > huge(x))  ! Does not always work
y = (.not. is_finite(x)) .and. (.not. is_inf(x))
end function is_nan_dp


#if __QP_AVAILABLE__ == 1

pure elemental function is_nan_qp(x) result(y)
use, non_intrinsic :: consts_mod, only : QP
implicit none
real(QP), intent(in) :: x
logical :: y
!y = (.not. (x <= huge(x) .and. x >= -huge(x))) .and. (.not. abs(x) > huge(x))  ! Does not always work
y = (.not. is_finite(x)) .and. (.not. is_inf(x))
end function is_nan_qp

#endif


end module infnan_mod
