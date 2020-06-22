      program testmod
      use consts, only : rp
      use lina
      implicit none
      
      real(kind = rp) :: A(2,3), B(3,2), C(3,1), D(1, 3), E(3), F(2, 1),&
     & G(1, 2), H(2), I(3), J(2)
      data A(1, :) /1, 2, 3/
      data A(2, :) /4, 5, 6/
      data B(1, :) /1, 2/
      data B(2, :) /3, 4/
      data B(3, :) /5, 6/
      data C(1, :) /1/
      data C(2, :) /2/
      data C(3, :) /3/ 
      data D(1, :) /1, 2, 3/
      E = [1, 2, 3]
      data F(1, :) /1/
      data F(2, :) /2/
      data G(1, :) /1, 2/
      H = [1, 2]
      I = [4, 5, 6]
      J = [3, 4]
!      print *, shape(matp(A,B)), size(matp(A,B))
!      print *, matp(A, B)
!      print *, shape(matp(A,C)), size(matp(A,C))
!      print *, matp(A, C)
      print *, matmul(A, B)
      print *, matmul(A, C)
      print *, matmul(G, A)
      print *, matmul(A, E)
      print *, matmul(H, A)
      print *, dot_product(E, I), dot_product(I, E)
      print *, dot_product(H, J), dot_product(J, H)
      print *, shape(outprod(E, I)), outprod(E, I)
      print *, shape(outprod(I,E)), outprod(I, E)
      print *, shape(outprod(H,E)), outprod(H, E)
      print *, shape(outprod(E,H)), outprod(E, H)
      end program testmod
