module matprod_mod

contains

function matprod22(x, y) result(z)

implicit none
integer, parameter :: RP = kind(0.0D0)
integer, parameter :: IK = kind(0)
real(RP), intent(in) :: x(:, :)
real(RP), intent(in) :: y(:, :)
real(RP) :: z(size(x, 1), size(y, 2))

integer(IK) :: i, j

z = 0.0D0
do j = 1, int(size(y, 2), kind(j))
    do i = 1, int(size(x, 2), kind(i))
        z(:, j) = z(:, j) + x(:, i) * y(i, j)
    end do
end do
end function matprod22
end module matprod_mod


program testmatmul

use matprod_mod, only : matprod22
implicit none
integer, parameter :: m = 2, n = 2, p = 2

real(kind(0.0D0)) :: A(m, n), B(n, p)

call random_number(A)
call random_number(B)
write (*, *) maxval(abs(matmul(A, B) - matprod22(A, B)))

end program testmatmul
