      module warnerror

      implicit none
      private
      public :: errmssg


      contains
          
      subroutine errmssg(srname, message)
      character(len=*), intent(in) :: srname
      character(len=*), intent(in) :: message 

      print *
      print '(1x, 5A)', 'Error: ',  srname(1 : len_trim(srname)),       &
     & ': ', message(1 : len_trim(message)), '.'
      print *

      end subroutine

      end module warnerror
