module ieee_4dev_mod
!--------------------------------------------------------------------------------------------------!
! This module makes some components of ieee_arithmetic available. Only for development and tests.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020
!
! Last Modified: Friday, March 04, 2022 PM05:23:55
!--------------------------------------------------------------------------------------------------!

use, intrinsic :: ieee_arithmetic, only : IEEE_VALUE, &
    & IEEE_QUIET_NAN, IEEE_SIGNALING_NAN, IEEE_POSITIVE_INF, IEEE_NEGATIVE_INF
use, non_intrinsic :: consts_mod, only : RP

implicit none

private
public :: ieeenan, ieeenan_q, ieeenan_s, ieeeinf, ieeeinf_p, ieeeinf_n


contains


pure real(RP) function ieeenan()
ieeenan = IEEE_VALUE(1.0_RP, IEEE_QUIET_NAN)
!ieeenan = IEEE_VALUE(1.0_RP, IEEE_SIGNALING_NAN)  ! Singling NaN can trigger an "floating invalid" error
end function ieeenan

pure real(RP) function ieeenan_q()
ieeenan_q = IEEE_VALUE(1.0_RP, IEEE_QUIET_NAN)
end function ieeenan_q

pure real(RP) function ieeenan_s()
ieeenan_s = IEEE_VALUE(1.0_RP, IEEE_SIGNALING_NAN)
end function ieeenan_s

pure real(RP) function ieeeinf()
ieeeinf = IEEE_VALUE(1.0_RP, IEEE_POSITIVE_INF)
end function ieeeinf

pure real(RP) function ieeeinf_p()
ieeeinf_p = IEEE_VALUE(1.0_RP, IEEE_POSITIVE_INF)
end function ieeeinf_p

pure real(RP) function ieeeinf_n()
ieeeinf_n = IEEE_VALUE(1.0_RP, IEEE_NEGATIVE_INF)
end function ieeeinf_n


end module ieee_4dev_mod
