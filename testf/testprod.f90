!        This is file : testprod
! Author= zaikunzhang
! Started at: 11.10.2021
! Last Modified: Monday, October 11, 2021 AM11:53:13

module test_mod
contains
function inprod(x, y) result(z)
real :: x(:), y(:)
real :: z
z = huge(0.0)
end function
end module test_mod

module testuse_mod
use test_mod, only : dot_product => inprod
use test_mod, only : inprod, inprod
end module testuse_mod

program testprod
use testuse_mod
implicit none

write (*, *) dot_product([1.0, 2.0], [3.0, 4.0])
write (*, *) inprod([1.0, 2.0], [3.0, 4.0])

end program testprod
