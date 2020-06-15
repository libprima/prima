      subroutine evalfun(n, x, f, xisnan)
      
      implicit none
      integer, parameter :: dp = kind(0.0d0)
      
      integer, intent(in) :: n
      real(kind = dp), intent(in) :: x(n)
      real(kind = dp), intent(out) :: f
      logical, intent(out) :: xisnan
      
      if (any(x /= x)) then
          xisnan = .true.
          f = sum(x)  ! set F to NaN; this is necessary 
      else
          xisnan = .false.
          call calfun(n, x, f)
      end if
      
      end subroutine evalfun
