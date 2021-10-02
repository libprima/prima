program testmaxmin
use, intrinsic :: iso_fortran_env, only : REAL128
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
implicit none

real(REAL128) :: a
write (*, *) minval([-1.0_REAL128, 1.0_REAL128, 0.0_REAL128])
write (*, *) minloc([-1.0_REAL128, 1.0_REAL128, 0.0_REAL128])
write (*, *) maxval([-1.0_REAL128, 1.0_REAL128, 0.0_REAL128])
write (*, *) maxloc([-1.0_REAL128, 1.0_REAL128, 0.0_REAL128])

a = 1.0E19_REAL128
write (*, *) minloc([a, a])
write (*, *) minloc([a, a], mask=[.true., .true.])

!write (*, *) ieee_is_nan(1.0_REAL128)

end program testmaxmin
