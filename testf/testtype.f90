!        This is file : testtype
! Author= zaikunzhang
! Started at: 08.11.2021
! Last Modified: Thursday, January 13, 2022 PM08:11:45

program testtype
use, intrinsic :: iso_fortran_env, only : INTEGER_KINDS
implicit none
integer, parameter :: ik = kind(1)
integer, parameter :: ikk = minval(INTEGER_KINDS, mask=INTEGER_KINDS /= ik)
integer(ik) :: i
integer(ikk) :: j

i = 0
j = 1

write (*, *) i, j


end program testtype
