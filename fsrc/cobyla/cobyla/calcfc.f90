subroutine calcfc(n, m, x, f, con)
! Test problem 10 (Hexagon area) in Powell's original COBYLA package.
!use, non_intrinsic :: consts_mod, only : ONE, HALF
!
implicit none

! Inputs
integer, parameter :: RP = kind(1.0)
integer :: m, n
real(RP), intent(in) :: x(n)


! Outputs
real(RP), intent(out) :: f
real(RP), intent(out) :: con(m)
real(RP) :: ONE = 1.0_RP
real(RP) :: HALF = 0.5_RP

f = -HALF * (x(1) * x(4) - x(2) * x(3) + x(3) * x(9) - x(5) * x(9) + x(5) * x(8) - x(6) * x(7))
con(1) = ONE - x(3)**2 - x(4)**2
con(2) = ONE - x(9)**2
con(3) = ONE - x(5)**2 - x(6)**2
con(4) = ONE - x(1)**2 - (x(2) - x(9))**2
con(5) = ONE - (x(1) - x(5))**2 - (x(2) - x(6))**2
con(6) = ONE - (x(1) - x(7))**2 - (x(2) - x(8))**2
con(7) = ONE - (x(3) - x(5))**2 - (x(4) - x(6))**2
con(8) = ONE - (x(3) - x(7))**2 - (x(4) - x(8))**2
con(9) = ONE - x(7)**2 - (x(8) - x(9))**2
con(10) = x(1) * x(4) - x(2) * x(3)
con(11) = x(3) * x(9)
con(12) = -x(5) * x(9)
con(13) = x(5) * x(8) - x(6) * x(7)
con(14) = x(9)

!write (*, *) x(1:n)
end subroutine calcfc
