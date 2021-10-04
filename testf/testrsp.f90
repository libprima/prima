!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! testrsp.f90

!!!!!! Module !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module test_mod
use, intrinsic :: iso_fortran_env, only : INT64
implicit none
private
public :: test

contains

subroutine test(ij)
implicit none
integer(INT64), intent(in) :: ij(:, :)

write (*, *) '*** In subroutine:'
write (*, *) 'Before reshaping: ', [maxval(ij, dim=2), minval(ij, dim=2)]
write (*, *) 'After reshaping: ', reshape([maxval(ij, dim=2), minval(ij, dim=2)], shape(ij))

end subroutine test

end module test_mod

!!!!!! Program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program testrsp
use, intrinsic :: iso_fortran_env, only : INT64
use, non_intrinsic :: test_mod, only : test
implicit none
integer(INT64) :: ij(1, 2)

ij(:, 1) = 1
ij(:, 2) = 2

! The following lines are OK. Benchmark for the subroutine, which is erroneous.
write (*, *) '*** In program:'
write (*, *) 'Before reshaping: ', [maxval(ij, dim=2), minval(ij, dim=2)]
write (*, *) 'After reshaping: ', reshape([maxval(ij, dim=2), minval(ij, dim=2)], shape(ij))

! The subroutine does exactly the same thing, but the result is wrong.
call test(ij)

end program testrsp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
