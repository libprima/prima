!        This is file : testqp
! Author= zaikunzhang
! Started at: 25.09.2021
! Last Modified: Saturday, September 25, 2021 PM09:54:25

program testqp
use, intrinsic :: iso_fortran_env, only : INT8, INT16, INT32, INT64, REAL32, REAL64, REAL128, INTEGER_KINDS, REAL_KINDS
implicit none

write (*, *) 'int8, int16, int32, int64', INT8, INT16, INT32, INT64

write (*, *) 'real32, real64, real128', REAL32, REAL64, REAL128

write (*, *) INTEGER_KINDS
write (*, *) REAL_KINDS


end program testqp
