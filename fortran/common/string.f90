module string_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some procedures for manipulating strings.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Monday, October 02, 2023 PM10:46:01
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: lower, upper, strip, istr, num2str

! MAX_NUM_STR_LEN is the maximum length of a string that is needed to represent a real or integer
! number. Assuming that such a number is represented by at most 128 bits, it is safe to set this
! maximum length to 128. We set this number to 256 to be on the safe side.
integer, parameter :: MAX_NUM_STR_LEN = 256
! MAX_WIDTH is the maximum number of characters printed in each row when printing arrays.
integer, parameter :: MAX_WIDTH = 100

interface num2str
    module procedure real2str_scalar, real2str_vector, int2str
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


function real2str_scalar(x, ndgt, nexp) result(s)
!--------------------------------------------------------------------------------------------------!
! This function converts a real scalar to a string. Optionally, NDGT is the number of decimal
! digits to print, and NEXP is the number of digits in the exponent; they may be reduced if needed.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DP, IK, REALMAX, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, validate
use, non_intrinsic :: infnan_mod, only : is_finite, is_nan
implicit none
! Inputs
real(RP), intent(in) :: x
integer, intent(in), optional :: ndgt
integer, intent(in), optional :: nexp
! Outputs
character(len=:), allocatable :: s
! Local variables
character(len=*), parameter :: srname = 'REAL2STR_SCALAR'
character(len=:), allocatable :: sformat
character(len=MAX_NUM_STR_LEN) :: str
integer(IK) :: ndgt_loc  ! The number of decimal digits to print
integer(IK) :: nexp_loc  ! The number of digits in the exponent
integer(IK) :: wx  ! The width of the printed X

! Preconditions
if (DEBUGGING) then
    if (present(ndgt)) then
        call assert(ndgt >= 0 .and. 2 * ndgt <= MAX_NUM_STR_LEN - 5, &
            & '0 <= NDGT <= '//int2str(floor(real(MAX_NUM_STR_LEN - 5) / 2.0, IK)), srname)
    end if
    if (present(nexp)) then
        call assert(nexp >= 0 .and. 2 * nexp <= MAX_NUM_STR_LEN - 5, &
            & '0 <= NEXP <= '//int2str(floor(real(MAX_NUM_STR_LEN - 5) / 2.0, IK)), srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

if (.not. is_finite(x)) then
    write (str, *) x
    s = strip(str)  ! Remove the leading and trailing spaces, if any.
else
    if (present(ndgt)) then
        ndgt_loc = int(ndgt, IK)
    else
        ! By default, we print at most the same number of decimal digits as the double precision.
        ndgt_loc = int(min(precision(x), precision(0.0_DP)), IK)
    end if
    ndgt_loc = min(ndgt_loc, floor(real(MAX_NUM_STR_LEN - 5) / 2.0, IK))
    if (present(nexp)) then
        nexp_loc = int(nexp, IK)
    else
        nexp_loc = ceiling(log10(real(range(x) + 0.1)), IK)  ! Use + 0.1 in case RANGE(X) = 10^k.
    end if
    nexp_loc = min(nexp_loc, floor(real(MAX_NUM_STR_LEN - 5) / 2.0, IK))
    wx = ndgt_loc + nexp_loc + 5_IK
    call validate(wx <= MAX_NUM_STR_LEN, 'The width of the printed number is at most ' &
        & //int2str(int(MAX_NUM_STR_LEN, IK)), srname)
    sformat = '(1PE'//int2str(wx)//'.'//int2str(ndgt_loc)//'E'//int2str(nexp_loc)//')'
    write (str, sformat) x
    s = trim(str)  ! Remove the trailing spaces, but keep the leading ones, if any.
end if

!====================!
!  Calculation ends  !
!====================!

! Postconditions
if (DEBUGGING) then
    call assert(len(s) > 0 .and. len(s) <= MAX_NUM_STR_LEN, '0 < LEN(S) <= MAX_NUM_STR_LEN', srname)
    call assert(is_nan(x) .eqv. is_nan(str2real(s)), 'IS_NAN(X) .EQV. IS_NAN(STR2REAL(S))', srname)
    ! The assertions concerning the infiniteness of X may fail due to the limited precision of
    ! printing. Thus we relax the assertions as below.
    !call assert(is_posinf(x) .eqv. is_posinf(str2real(s)), 'IS_POSINF(X) .EQV. IS_POSINF(STR2REAL(S))', srname)
    !call assert(is_neginf(x) .eqv. is_neginf(str2real(s)), 'IS_NEGINF(X) .EQV. IS_NEGINF(STR2REAL(S))', srname)
    call assert(x >= REALMAX * (1.0 - 10.0**(-ndgt_loc)) .eqv. str2real(s) >= REALMAX * (1.0 - 10.0**(-ndgt_loc)), &
        & 'IS_POSINF(X) .EQV. IS_POSINF(STR2REAL(S))', srname)
    call assert(x <= -REALMAX * (1.0 - 10.0**(-ndgt_loc)) .eqv. str2real(s) <= -REALMAX * (1.0 - 10.0**(-ndgt_loc)), &
        & 'IS_NEGINF(X) .EQV. IS_NEGINF(STR2REAL(S))', srname)
    if (abs(x) < REALMAX) then
        call assert(abs(x - str2real(s)) <= abs(x) * 10.0**(-ndgt_loc), 'STR2REAL(S) == X', srname)
    end if
end if
end function real2str_scalar

function real2str_vector(x, ndgt, nexp, nx) result(s)
!--------------------------------------------------------------------------------------------------!
! This function converts a real vector to a string. Optionally, NDGT is the number of decimal
! digits to print, NEXP is the number of digits in the exponent, and NX is the number of entries
! printed per row; they may be reduced if needed.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: memory_mod, only : safealloc
! Inputs
real(RP), intent(in) :: x(:)
integer, intent(in), optional :: ndgt
integer, intent(in), optional :: nexp
integer, intent(in), optional :: nx
! Outputs
character(len=:), allocatable :: s
! Local variables
character(len=*), parameter :: srname = 'REAL2STR_VECTOR'
character(len=2), parameter :: spaces = '  '  ! The spaces between two entries in a row
integer :: i
integer :: j
integer :: m  ! The number of rows
integer :: n  ! N = SIZE(X)
integer :: ndgt_loc  ! The number of decimal digits to print
integer :: nexp_loc  ! The number of digits in the exponent
integer :: nx_loc  ! The number of entries printed per row
integer :: slen  ! The length of the string
integer :: wx  ! The width of each entry in X

! Preconditions
if (DEBUGGING) then
    if (present(ndgt)) then
        call assert(ndgt >= 0 .and. 2 * ndgt <= MAX_NUM_STR_LEN - 5, &
            & '0 <= NDGT <= '//int2str(floor(real(MAX_NUM_STR_LEN - 5) / 2.0, IK)), srname)
    end if
    if (present(nexp)) then
        call assert(nexp >= 0 .and. 2 * nexp <= MAX_NUM_STR_LEN - 5, &
            & '0 <= NEXP <= '//int2str(floor(real(MAX_NUM_STR_LEN - 5) / 2.0, IK)), srname)
    end if
    if (present(nx)) then
        call assert(nx >= 1, 'NX >= 1', srname)
    end if
end if

!====================!
! Calculation starts !
!====================!

! Quick return if X is empty.
if (size(x) <= 0) then
    s = ''
    return
end if

if (present(ndgt)) then
    ndgt_loc = ndgt
else
    ! By default, we print at most the same number of decimal digits as the double precision.
    ndgt_loc = min(precision(x), precision(0.0_DP))
end if
ndgt_loc = min(ndgt_loc, floor(real(MAX_NUM_STR_LEN - 5) / 2.0))

if (present(nexp)) then
    nexp_loc = nexp
else
    nexp_loc = ceiling(log10(real(range(x)) + 0.1))  ! Use + 0.1 in case RANGE(X) = 10^k.
end if
nexp_loc = min(nexp_loc, floor(real(MAX_NUM_STR_LEN - 5) / 2.0))

wx = len(real2str_scalar(0.0_RP, ndgt, nexp))
n = size(x)
if (present(nx)) then
    nx_loc = max(1, min(nx, n))
else
    nx_loc = max(1, min(floor(real(MAX_WIDTH + len(spaces)) / (real(wx) + len(spaces))), size(x)))
end if

! Calculate the length of the printed string S.
! N.B.: Here, SLEN should not be an INT16 integer, because a double-precision vector of length
! ~3500 would be printed as a string longer than 65536. On most modern platforms, the default
! integer kind is INT32, which is enough for printing double-precision vectors of size ~ 10^8,
! being sufficient for this project.
m = ceiling(real(n) / real(nx_loc))  ! The number of rows
slen = wx * n + len(spaces) * (n - 1) + (1 - len(spaces)) * (m - 1)
call safealloc(s, slen)

j = 0  ! J is the index of the last up-to-date character in S.
do i = 1, n
    s(j + 1:j + wx) = real2str_scalar(x(i), ndgt_loc, nexp_loc)
    if (i == n) exit
    j = j + wx
    if (modulo(i, nx_loc) == 0) then
        s(j + 1:j + 1) = new_line(s)
        j = j + 1
    else
        s(j + 1:j + len(spaces)) = spaces
        j = j + len(spaces)
    end if
end do

!====================!
!  Calculation ends  !
!====================!
end function real2str_vector

function str2real(s) result(x)
!--------------------------------------------------------------------------------------------------!
! This function converts a string to a real scalar.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none
character(len=*), intent(in) :: s
real(RP) :: x
character(len=*), parameter :: srname = 'STR2REAL'
if (DEBUGGING) then
    call assert(len(s) > 0, 'LEN(S) > 0', srname)
end if
read (s, *) x
end function str2real


function int2str(x) result(s)
!--------------------------------------------------------------------------------------------------!
! This function converts an integer scalar to a string.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none
integer(IK), intent(in) :: x
character(len=*), parameter :: srname = 'INT2STR'
character(len=:), allocatable :: s
character(len=MAX_NUM_STR_LEN) :: str
write (str, *) x
s = strip(str)
if (DEBUGGING) then
    call assert(len(s) > 0 .and. len(s) <= MAX_NUM_STR_LEN, '0 < LEN(S) <= MAX_NUM_STR_LEN', srname)
    call assert(str2int(s) == x, 'STR2INT(S) == X', srname)
end if
end function int2str

function str2int(s) result(x)
!--------------------------------------------------------------------------------------------------!
! This function converts a string to an integer scalar.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
implicit none
character(len=*), intent(in) :: s
integer(IK) :: x
character(len=*), parameter :: srname = 'STR2INT'
if (DEBUGGING) then
    call assert(len(s) > 0, 'LEN(S) > 0', srname)
end if
read (s, *) x
end function str2int


end module string_mod
