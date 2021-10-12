program testk

use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_support_nan
implicit none

if (.not. ieee_support_nan(1.)) error stop "Lack of support"
print *, kind(ieee_is_nan(1.)), kind(.true.)

end program testk
