program testov
use, intrinsic :: iso_fortran_env, only : REAL128
implicit none

real(REAL128) :: a
logical :: b

a = -huge(0.0_REAL128)
b = (a < huge(0.0_REAL128))

write (*, *) b

end program testov
