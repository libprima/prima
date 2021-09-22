!        This is file : testgeneric
! Author= zaikunzhang
! Started at: 14.09.2021
! Last Modified: Tuesday, September 14, 2021 AM10:54:44

module generic_mod

implicit none
private
public :: string

contains

pure function string(numeric) result(number_as_string)
class(*), intent(in) :: numeric
integer, parameter :: max_len = 128
character(len=max_len) :: untrimmed_string
character(len=:), allocatable :: number_as_string

select type (numeric)
type is (complex)
    write (untrimmed_string, *) numeric
type is (integer)
    write (untrimmed_string, *) numeric
type is (logical)
    write (untrimmed_string, *) numeric
type is (real)
    write (untrimmed_string, *) numeric
type is (real(kind(0.0D0)))
    write (untrimmed_string, *) numeric
type is (real(-1))
class default
    error stop "Internal error in subroutine 'assert': unsupported type in function 'string'."
end select


number_as_string = trim(adjustl(untrimmed_string))

end function string

end module generic_mod


program testgeneric
use generic_mod, only : string
implicit none
write (*, *) string(1)
write (*, *) string(1.0)
write (*, *) string(1.0D0)
write (*, *) string(.true.)

end program testgeneric
