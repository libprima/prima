!        This is file : test_allocstr
! Author= zaikunzhang
! Started at: 05.12.2021
! Last Modified: Sunday, December 05, 2021 PM04:09:07

program test_allocstr
implicit none
character(:), allocatable :: a
allocate (a(2))
a = 'ab'
write (*, *) a

end program test_allocstr
