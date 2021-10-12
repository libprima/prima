program testany
implicit none
logical :: A(2, 3)
A(1, :) = [.true., .true., .false.]
A(2, :) = [.false., .true., .true.]
write (*, *) all(A)
end program testany
