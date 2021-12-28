module datetime_mod
implicit none
private
public :: year
public :: week


contains


function year(ymd) result(y)
!--------------------------------------------------------------------------------------------------!
! This function returns the number of the current year.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: ymd(3)
! Outputs
integer :: y

! Local values
integer :: values(8)

if (present(ymd)) then
    y = ymd(1)
else
    call date_and_time(values=values)
    y = values(1)
end if
end function year


!function month(ymd) result(m)
!!--------------------------------------------------------------------------------------------------!
!! This function returns the number of the current month.
!!--------------------------------------------------------------------------------------------------!
!implicit none
!! Inputs
!integer, intent(in), optional :: ymd(3)
!! Outputs
!integer :: m

!! Local values
!integer :: values(8)

!if (present(ymd)) then
!    m = ymd(2)
!else
!    call date_and_time(values=values)
!    m = values(2)
!end if
!end function month


function day(ymd) result(d)
!--------------------------------------------------------------------------------------------------!
! This function returns the number of the day YMD(3) in the year YMD(1) or the current day in the
! current year.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: ymd(3)
! Outputs
integer :: d

! Local values
integer :: dom
integer :: monlens(12)
integer :: moy
integer :: values(8)
integer :: y

if (present(ymd)) then
    y = ymd(1)
    moy = ymd(2)
    dom = ymd(3)
else
    call date_and_time(values=values)
    y = values(1)
    moy = values(2)
    dom = values(3)
end if

monlens = month_lengths(y)
d = sum(monlens(1:moy - 1)) + dom
end function day


function week(ymd) result(w)
!--------------------------------------------------------------------------------------------------!
! This function returns the number of the current week in the current year.
! N.B.: This function takes the week of January 1 as the first week and increments W every MONDAY.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: ymd(3)
! Outputs
integer :: w

! Local variables
integer :: jan1
integer :: values(8)
integer :: y

if (present(ymd)) then
    y = ymd(1)
else
    call date_and_time(values=values)
    y = values(1)
end if

jan1 = whatday([y, 1, 1])

if (present(ymd)) then
    w = ceiling(real(day(ymd) + jan1 - 1) / 7.0)
else
    w = ceiling(real(day() + jan1 - 1) / 7.0)
end if
end function week


function isleap(y) result(is_leap)
!--------------------------------------------------------------------------------------------------!
! This function tests whether the year Y or the current year is a leap year.
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


function days_since_20000101(ymd) result(d)
!--------------------------------------------------------------------------------------------------!
! This function returns the number of days since 2000-01-01 until YMD or the current day.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: ymd(3)
! Outputs
integer :: d

integer :: values(12)
integer :: y

if (present(ymd)) then
    y = ymd(1)
else
    call date_and_time(values=values)
    y = values(1)
end if

d = 0
if (y > 2000) then
    do y = 2000, y - 1
        if (isleap(y)) then
            d = d + 366
        else
            d = d + 365
        end if
    end do
elseif (y < 2000) then
    do y = y, 1999
        if (isleap(y)) then
            d = d - 366
        else
            d = d - 365
        end if
    end do
end if
d = d + day(ymd) - 1
end function days_since_20000101


function whatday(ymd) result(d)
!--------------------------------------------------------------------------------------------------!
! This function returns the day of week of the day YMD(3) or the current day, 1 for Monday, and
! 7 for Sunday.
! N.B.: 2000-01-01 is Saturday.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: ymd(3)
! Outputs
integer :: d

if (present(ymd)) then
    d = mod(5 + days_since_20000101(ymd), 7) + 1
else
    d = mod(5 + days_since_20000101(), 7) + 1
end if
if (d <= 0) then
    d = d + 7
end if
end function whatday


function month_lengths(y) result(monlens)
!--------------------------------------------------------------------------------------------------!
! This function returns the lengths of months of the year Y or the current year.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
integer, intent(in), optional :: y
! Outputs
integer :: monlens(12)

if (isleap(y)) then
    monlens = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
else
    monlens = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
end if
end function month_lengths


end module datetime_mod
