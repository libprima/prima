program test

use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: prob_mod, only : getdelta0

implicit none
integer(IK) :: n

n = 1

write (*, *) getdelta0(n)

end program test
