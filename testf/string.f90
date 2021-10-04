module string_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some procedures for manipulating strings.
!
! Coded by Zaikun Zhang in September 2021.
!
! Last Modified: Monday, September 20, 2021 PM07:53:25
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: lower, upper, trimstr


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


end module string_mod
