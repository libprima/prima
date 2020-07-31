
        module overload

        interface sqrt
            module procedure char_sqrt
        end interface

        contains

        function char_sqrt (x)
        real :: char_sqrt
        character(*), intent(in) :: x
        real :: xr
        read (x,*) xr
        char_sqrt = sqrt(xr)
        end function char_sqrt
        end module overload

        program test
        use overload, only : sqrt
        print *, sqrt('2.0'), sqrt(2.0)
        end program test 
