module datetime_mod
implicit none
private
public :: year
public :: week
public :: isleap


contains


function year() result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns the number of the current year.
!--------------------------------------------------------------------------------------------------!
implicit none
! Outputs
integer :: y

! Local values
integer :: values(8)

call date_and_time(values=values)
y = values(1)
end function year


!function month() result(m)
!--------------------------------------------------------------------------------------------------!
! This function returns the number of the current month.
!--------------------------------------------------------------------------------------------------!
!implicit none
!! Outputs
!integer :: m

!! Local values
!integer :: values(8)

!call date_and_time(values=values)
!m = values(2)
!end function month


function day(ymd) result(d)
!--------------------------------------------------------------------------------------------------!
! This function returns the number of the current day in the current year.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: ymd(3)
! Outputs
integer :: d

! Local values
integer :: dom
integer :: monthday(12)
integer :: moy
integer :: values(8)
integer :: year

if (present(ymd)) then
    year = ymd(1)
    moy = ymd(2)
    dom = ymd(3)
else
    call date_and_time(values=values)
    year = values(1)
    moy = values(2)
    dom = values(3)
end if

if (isleap(year)) then
    monthday = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
else
    monthday = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
end if
d = sum(monthday(1:moy - 1)) + dom
end function day


function week(ymd) result(w)
!--------------------------------------------------------------------------------------------------!
! This function returns the number of the current week in the current year.
! N.B.: This function increments W every 7 days starting from January 1, not necessarily on Monday.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: ymd(3)
! Outputs
integer :: w

if (present(ymd)) then
    w = ceiling(real(day(ymd)) / 7.0)
else
    w = ceiling(real(day()) / 7.0)
end if
end function week


function isleap(y) result(is_leap)
!--------------------------------------------------------------------------------------------------!
! This function tests whether the current year is a leap year.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: y
! Outputs
logical :: is_leap

! Local variables
integer :: y_loc

if (present(y)) then
    y_loc = y
else
    y_loc = year()
end if

is_leap = (mod(y_loc, 4) == 0 .and. mod(y_loc, 100) /= 0) .or. (mod(y_loc, 400) == 0)
end function isleap

end module datetime_mod
