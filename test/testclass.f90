program testclass
write (*, *) string(2)
end program testclass

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
class default
    error stop "Internal error in subroutine 'assert': unsupported type in function 'string'."
end select

number_as_string = trim(adjustl(untrimmed_string))

end function string
