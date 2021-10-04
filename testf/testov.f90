program testov
use, intrinsic :: iso_fortran_env, only : REAL128
use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
implicit none

real(REAL128) :: a
logical :: b

a = -huge(0.0_REAL128)
b = (a < huge(0.0_REAL128))

write (*, *) b
!write (*, *) ieee_is_nan(a)
write (*, *) (.not. (a <= huge(a) .and. a >= -huge(a))) .and. (.not. abs(a) > huge(a))
!write(*,*) a>= huge(a)
!write(*,*) a<= huge(a)

write (*, *) minval([-1.0_REAL128, 1.0_REAL128, 0.0_REAL128])
write (*, *) minloc([-1.0_REAL128, 1.0_REAL128, 0.0_REAL128])
write (*, *) maxval([-1.0_REAL128, 1.0_REAL128, 0.0_REAL128])
write (*, *) maxloc([-1.0_REAL128, 1.0_REAL128, 0.0_REAL128])


end program testov
