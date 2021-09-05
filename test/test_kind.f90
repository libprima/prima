! test_kind.f90
program test_kind

use ieee_arithmetic, only : ieee_is_nan
implicit none

write (*, *) kind(1)
write (*, *) kind(ieee_is_nan(1.0))  ! This prints 4.
write (*, *) kind(.false.)   ! A comparison, which prints 4.

end program test_kind
