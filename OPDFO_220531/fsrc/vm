module vm_mod

implicit none
private
public :: v2m, m2v


contains


function v2m(hqv) result(hqm)
use, non_intrinsic :: consts_mod, only : RP, IK, ZERO
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : issymmetric
implicit none
real(RP), intent(in) :: hqv(:)
real(RP) :: hqm((floor(sqrt(real(8 * size(hqv) + 1))) - 1) / 2, (floor(sqrt(real(8 * size(hqv) + 1))) - 1) / 2)
integer(IK) :: n, i, j, ih

n = int(size(hqm, 1), kind(n))
call assert(size(hqv) == n * (n + 1) / 2, 'SIZE(HQV) = N*(N+1)/2', 'v2m')

hqm = ZERO
ih = 0_IK
do j = 1_IK, n
    do i = 1_IK, j
        ih = ih + 1_IK
        hqm(i, j) = hqv(ih)
        hqm(j, i) = hqv(ih)
    end do
end do

call assert(issymmetric(hqm), 'HQM is symmetric', 'V2M')

end function v2m

function m2v(hqm) result(hqv)
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : issymmetric
implicit none
real(RP), intent(in) :: hqm(:, :)
real(RP) :: hqv((size(hqm, 1) * (size(hqm, 1) + 1)) / 2)
integer(IK) :: n, i, j, ih

call assert(issymmetric(hqm), 'HQM is symmetric', 'M2V')

n = int(size(hqm, 1), kind(n))
ih = 0_IK
do j = 1_IK, n
    do i = 1_IK, j
        ih = ih + 1_IK
        hqv(ih) = hqm(i, j)
    end do
end do

end function m2v

end module vm_mod
