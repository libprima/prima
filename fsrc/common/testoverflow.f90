!        This is file : testoverflow
! Author= zaikunzhang
! Started at: 16.11.2021
! Last Modified: Tuesday, November 16, 2021 AM11:51:07

program testoverflow
implicit none
integer, parameter :: RK = kind(0.0)
real, parameter :: x = sqrt(huge(0.0))
real, parameter :: y = x * 3.141592653589793238462643383279502884_RK
real, parameter :: z(2) = [y, y]

end program testoverflow
