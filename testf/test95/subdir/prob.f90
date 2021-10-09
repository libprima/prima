module prob_mod
implicit none
private
public :: getdelta0

contains


function getdelta0(n) result(delta)

use, non_intrinsic :: consts_mod, only : IK

implicit none

integer(IK), intent(in) :: n
real :: delta

delta = 1.0


end function getdelta0

end module prob_mod
