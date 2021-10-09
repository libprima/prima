program testnan

use, intrinsic :: iso_fortran_env, only : QP => REAL128

implicit none

real(QP) :: x = 52936146080073098.443184315494389_QP

write (*, *) cos(x)
write (*, *) cos(cos(x))

end program testnan
