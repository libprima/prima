module kindt_mod
      implicit none
      integer, parameter :: IK = kind(0)

      interface 
          subroutine kindt_sub(i)
          implicit none
          integer(IK), intent(in) :: i
          end subroutine kindt_sub
      end interface
!      contains
!      subroutine kindt_sub(i)
!      integer(IK), intent(in) :: i
!      print *, i
!      end subroutine kindt_sub
end module kindt_mod

program kindt
      implicit none
      use kindt_mod, only : kindt_sub
      call kindt_sub(1)
end program kindt

subroutine kindt_sub(i)
implicit none
use kindt_sub, only : IK
integer(IK), intent(in) :: i
print *, i
end subroutine kindt_sub
