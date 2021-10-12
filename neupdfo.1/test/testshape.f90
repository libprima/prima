program testshape
real, external :: test(:, :)
real :: A(3, 4)
real :: B(4, 3)
B = test(A)
end program testshape


function test(A) result(B)
real :: A(:, :)
real :: B(size(A, 2), size(A, 1))
B = transpose(A)
end function test
