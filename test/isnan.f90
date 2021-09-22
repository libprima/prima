module isnan_mod

implicit none
private
public :: is_nan

contains

pure function is_nan(x) result(y)

class(*), intent(in) :: x

logical :: y

select type (x)
type is (real)
    y = (.not. x > 0) .and. (.not. x < 1)
type is (real(kind(0.0D0)))
    y = (.not. x > 0) .and. (.not. x < 1)
class default
    y = .true.
end select

end function is_nan


end module isnan_mod
