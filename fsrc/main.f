      program testmod
      use consts, only : rp
      use lina
      implicit none
      
      integer :: ind(2), jnd(2)
      real(kind = rp) :: A(2,3), B(3,2), C(3,1), D(1, 3), E(3), F(2, 1),&
     & G(1, 2), H(2), I(3), J(2), s, co, r
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
      ind = [1,3]
      jnd = [1,2]
      !print *, shape(A(jnd, ind)), A(jnd, ind) 
      call givens(0.0d0, 0.0d0, s, co, r)
      print *, s, co, r
      call givens(1.0d0, 0.0d0, s, co, r)
      print *, s, co, r
      call givens(-1.0d0, 0.0d0, s, co, r)
      print *, s, co, r
      call givens(0.0d0, -1.0d0, s, co, r)
      print *, s, co, r
      call givens(0.0d0, 1.0d0, s, co, r)
      print *, s, co, r
      call givens(1.0d0, 1.0d0, s, co, r)
      print *, s, co, r
      call givens(-1.0d0, 2.0d0, s, co, r)
      print *, s, co, r
      call givens(-1.0d0, -2.0d0, s, co, r)
      print *, s, co, r
      call givens(1.0d0, -2.0d0, s, co, r)
      print *, s, co, r
      call givens(1.0d0, 2.0d0, s, co, r)
      print *, s, co, r
      call givens(2.0d0, 1.0d0, s, co, r)
      print *, s, co, r
      call givens(-2.0d0, -1.0d0, s, co, r)
      print *, s, co, r
      call givens(2.0d0, -1.0d0, s, co, r)
      print *, s, co, r
      end program testmod
