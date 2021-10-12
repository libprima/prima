program testmem
      use iso_fortran_env, only : INT64
      integer(INT64) :: l 
      real(kind(0.0D0)), allocatable :: x(:)
      l = 2_INT64**35 
      l = l/8
      do i = 1, 100
          allocate(x(l))
          deallocate(x)
          print *, i
      end do
end program testmem
!program test
!   implicit none
!   real,allocatable :: array(:) 
!
!   allocate(array(1000000000)) !4 gb array
!
!   print *,'Check memory - 4 GB allocated'
!   read *
!
!   array(1:1000000) = 1.0
!
!   print *,'Check memory - 4 MB assigned'
!   read *
!
!   array(1000000:100000000) = 2.0
!
!   print *,'Check memory - 400 MB assigned'
!   read *
!
!   array = 5.0
!
!   print *,'Check memory - 4 GB assigned'
!   read *
!
!end program
