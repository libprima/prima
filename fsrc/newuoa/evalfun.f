      subroutine evalfun(n, x, f, xisnan)
      
      use pdfomod, only : rp, is_nan
      implicit none
      
      integer, intent(in) :: n
      real(kind = rp), intent(in) :: x(n)
      real(kind = rp), intent(out) :: f
      logical, intent(out) :: xisnan
      
      if (any(is_nan(x))) then
          xisnan = .true.
          f = sum(x)  ! set F to NaN; this is necessary 
      else
          xisnan = .false.
          call calfun(n, x, f)
      end if
      
      end subroutine evalfun
