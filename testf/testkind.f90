!        This is file : testkind
! Author= zaikunzhang
! Started at: 09.10.2021
! Last Modified: Saturday, October 09, 2021 PM05:30:11
!
module test_mod
implicit none
integer, parameter :: INT16 = selected_int_kind(3)
integer, parameter :: INT32 = selected_int_kind(7)
integer, parameter :: INT64 = selected_int_kind(14)
integer, parameter :: INT_DFT = kind(0)
contains
function test(n) result(y)
integer(INT16), intent(in) :: n
integer(INT16) :: y
y = n
end function test
end module test_mod

program testkind
!use, intrinsic :: iso_fortran_env, only : INT16, INT32, INT64
use test_mod, only : test
implicit none

integer, parameter :: INT16 = selected_int_kind(3)
integer, parameter :: INT32 = selected_int_kind(7)
integer, parameter :: INT64 = selected_int_kind(14)
integer, parameter :: INT_DFT = kind(0)

write (*, *) INT16, INT32, INT64, INT_DFT
write (*, *) huge(0_INT16), huge(0_INT32), huge(0_INT64), huge(0_INT_DFT), huge(0)

write (*, *) test(0_INT16)
end program testkind
