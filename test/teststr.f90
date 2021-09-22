! teststr.f90

!!!!!!! A module that defines COPY !!!!!!!!!!!!!!!!!!!!!!!!
module string_mod

implicit none
private
public :: copy

contains

pure function copy(x) result(y)
! A function that does nothing but copying X to Y.

implicit none
character(len=*), intent(in) :: x
character(len=len(x)) :: y
y = x
end function copy

end module string_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!! The main program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program teststr
! A simple test of COPY.
use, non_intrinsic :: string_mod, only : copy
implicit none

write (*, *) copy('a')//'b'

! The following line will not cause the warning.
!!write (*, *) copy('abc')
end program teststr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
