! testshape.f90
program testshape
implicit none
integer :: a(3, 2)
a(:, :) = a(:, :)  ! This line is OK. Benchmark for the erroneous line.
a(:, 1:2) = a(:, 1:2)  ! This line is OK. Benchmark for the erroneous line.
a(:, [1, 2]) = a(:, [1, 2])  ! Do nothing ..., but it triggers the error.
end program testshape
