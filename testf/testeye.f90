module eye_mod
use, intrinsic :: iso_fortran_env, only : INTEGER_KINDS
implicit none

integer, parameter :: IK = INTEGER_KINDS(1) ! or whatever kind that is valid
integer, parameter, private :: ID = kind(1)
logical, parameter, private :: tf = (INTEGER_KINDS(1) == IK)
integer, private :: ii = merge(2, 1, tf)
integer, parameter, private :: IOTHER = merge(ID, INTEGER_KINDS(ii), ID /= IK)

!interface eye
!    module procedure eye_ik, eye_iother
!end interface eye

!contains

!function eye_ik(n)
!implicit none
!integer(IK), intent(in) :: n
!real :: eye_ik(n, n)
!integer(IK) :: j, k
!eye_ik = reshape([1.0, ([(0.0, k=1, n)], 1.0, j=1, n - 1)], shape(eye_ik))
!end function

!function eye_iother(n)
!implicit none
!integer(IOTHER), intent(in) :: n
!real :: eye_iother(n, n)
!integer(IOTHER) :: j, k
!eye_iother = reshape([1.0, ([(0.0, k=1, n)], 1.0, j=1, n - 1)], shape(eye_iother))
!end function

end module eye_mod

!program testeye
!use eye_mod
!implicit none
!print *, eye(3)
!print *, eye(3_IK)
!end program testeye
