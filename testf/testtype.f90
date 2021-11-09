!        This is file : testtype
! Author= zaikunzhang
! Started at: 08.11.2021
! Last Modified: Tuesday, November 09, 2021 AM10:18:58

program testtype
use, intrinsic :: iso_fortran_env, only : INTEGER_KINDS
implicit none
integer, parameter :: ik = kind(1)
integer, parameter :: ikk = minval(INTEGER_KINDS, mask=INTEGER_KINDS /= ik)
integer, parameter :: ikkk = 21
integer(ik) :: i
integer(ikk) :: j
integer(ikkk) :: k

i = 0
j = 1
k = 2

write (*, *) i, j, k


end program testtype
