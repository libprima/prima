module string_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some procedures for manipulating strings.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Monday, April 10, 2023 PM12:19:56
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: lower, upper, strip, istr, num2str

! MAX_NUM_STR_LEN is the maximum length of a string that is needed to represent a real or integer
! number. Assuming that such a number is represented by at most 128 bits, it is safe to set this
! maximum length to 128.
integer, parameter :: MAX_NUM_STR_LEN = 128

interface num2str
    module procedure real2str, int2str
end interface num2str


contains


pure function lower(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function maps the characters of a string to the lower case, if applicable.
!--------------------------------------------------------------------------------------------------!

implicit none

character(len=*), intent(in) :: x
character(len=len(x)) :: y

integer, parameter :: dist = ichar('A') - ichar('a')
integer :: i

y = x
do i = 1, len(y)
    if (y(i:i) >= 'A' .and. y(i:i) <= 'Z') then
        y(i:i) = char(ichar(y(i:i)) - dist)
    end if
end do
end function lower


pure function upper(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function maps the characters of a string to the upper case, if applicable.
!--------------------------------------------------------------------------------------------------!

implicit none

character(len=*), intent(in) :: x
character(len=len(x)) :: y

integer, parameter :: dist = ichar('A') - ichar('a')
integer :: i

y = x
do i = 1, len(y)
    if (y(i:i) >= 'a' .and. y(i:i) <= 'z') then
        y(i:i) = char(ichar(y(i:i)) + dist)
    end if
end do
end function upper


pure function strip(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function removes the leading and trailing spaces of a string.
!--------------------------------------------------------------------------------------------------!

implicit none

character(len=*), intent(in) :: x
character(len=len(trim(adjustl(x)))) :: y

y = trim(adjustl(x))
end function strip


pure function istr(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function converts a string to an integer array.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
implicit none
character(len=*), intent(in) :: x
integer(IK) :: y(len(x))

integer(IK) :: i

y = [(int(ichar(x(i:i)), IK), i=1, int(len(x), IK))]

end function istr


pure function real2str(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function converts a real number to a string.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP
implicit none
real(RP), intent(in) :: x
character(len=:), allocatable :: y
character(len=MAX_NUM_STR_LEN) :: str
write (str, *) x
y = strip(str)
end function real2str


pure function int2str(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function converts an integer number to a string.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
implicit none
integer(IK), intent(in) :: x
character(len=:), allocatable :: y
character(len=MAX_NUM_STR_LEN) :: str
write (str, *) x
y = strip(str)
end function int2str


end module string_mod
