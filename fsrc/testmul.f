      program testmul
          implicit none
          real :: A(3,4), B(3,2)
          integer :: n

          A(1, :) = [1, 2, 3, 4]
          A(2, :) = [5, 6, 7, 8]
          A(3, :) = [9, 10, 11, 12]
          B(1, :) = [1, 2]
          B(2, :) = [3, 4]
          B(3, :) = [5, 6]

          print *, matmul(A(:, 1:3), B)

          call setn(n)

          print *, n
          print *, matmul(A(:, 1:n), B)
          print *, size(matmul(A(:, 1:n), B), 1)
          print *, size(matmul(A(:, 1:n), B), 2)

          contains 
          subroutine setn(n)
            implicit none
            integer ::n
            n = 5 
          end subroutine


      end program testmul

      module mul
          
      function matmul22(x, y) result(z)
      implicit none
      real(kind = rp), intent(in) :: x, y
      real(kind = rp), dimension(size(x, 1), size(y, 2)) :: z
      integer :: i, j
      if (size(x, 2) /= size(y, 1)) then
          stop 'Error in MATMUL: SIZE(X, 2) /= SIZE(Y, 1).'
      end if
      z = zero
      do j = 1, size(y, 2)
          do i = 1, size(x, 2)
              z(:, j) = z(:, j) + x(:, i)*y(i, j)
          end do
      end do
      end function matmul22
      end module mul
