module string_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some procedures for manipulating strings.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Wednesday, September 22, 2021 PM06:52:11
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: lower, upper, trimstr, istr


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


pure function trimstr(x) result(y)
!--------------------------------------------------------------------------------------------------!
! This function removes the leading and trailing spaces of a string.
!--------------------------------------------------------------------------------------------------!

implicit none

character(len=*), intent(in) :: x
character(len=len(trim(adjustl(x)))) :: y

y = trim(adjustl(x))
end function trimstr


pure function istr(x) result(y)
use, non_intrinsic :: consts_mod, only : IK
implicit none
character(len=*), intent(in) :: x
integer(IK) :: y(len(x))

integer(IK) :: i

y = [(int(ichar(x(i:i)), IK), i=1, int(len(x), IK))]

end function istr

end module string_mod
