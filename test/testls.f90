!        This is file : testls
! Author= zaikunzhang
! Started at: 22.09.2021
! Last Modified: Wednesday, September 22, 2021 AM10:31:41

module kind_mod
implicit none
private
public :: DP
integer, parameter :: DP = kind(0.0D0)
end module kind_mod

program testls
use kind_mod, only : DP
implicit none
real(DP) :: x
x = 0.0_DP
end program testls
