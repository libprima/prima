!        This is file : testcos
! Author= zaikunzhang
! Started at: 09.09.2021
! Last Modified: Thursday, September 09, 2021 AM11:01:02

program testcos
use iso_fortran_env, only : REAL64, REAL128
implicit none
write (*, *) cos(0.59843577329095299_REAL64)
write (*, *) cos(0.59843577329095299_REAL128)

end program testcos
