      
      program test

      use act_mod
      implicit none

      interface 
        subroutine calfun(f, x, y)

      integer, intent(in) :: x
      integer, intent(in) ::y 
      integer, intent(out) :: f
        end subroutine

      end interface 

       

      call act(calfun, 10)
      

      end program

      subroutine calfun(f, x, y)
      integer, intent(in) :: x
      integer, intent(in) ::y 
      integer, intent(out) :: f

      f = x-y

      end subroutine
