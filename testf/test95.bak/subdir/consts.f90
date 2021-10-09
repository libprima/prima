module consts_mod

implicit none
private
public :: IK

integer, parameter :: INT32 = selected_int_kind(7)

integer, parameter :: IK = INT32

end module consts_mod
