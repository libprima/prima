      subroutine update(n, npt, bmat, zmat, idz, vlag, beta, knew)
      ! UPDATE updates arrays BMAT and ZMAT together with IDZ, in order
      ! to shift the interpolation point that has index KNEW. On entry, 
      ! VLAG contains the components of the vector THETA*WCHECK + e_b 
      ! of the updating formula (6.11) in the NEWUOA paper, and BETA
      ! holds the value of the parameter that has this name. 

      use consts, only : one, zero
      use lina
      implicit none

      integer, intent(in) :: n, npt, knew
      integer, intent(inout) :: idz
      real(kind=rp), intent(in) :: beta
      real(kind=rp), intent(inout) :: bmat(npt+n, n), zmat(npt,npt-n-1),&
     & vlag(npt + n)

      integer :: iflag, j, ja, jb, jl
      real(kind=rp) :: c, s, r, alpha, denom, scala, scalb, tau, tausq, &
     & temp, tempa, tempb, ztemp(npt), w(npt + n), v1(n), v2(n)

      ! Apply the rotations that put zeros in the KNEW-th row of ZMAT.
      do j = 1, npt - n - 1
          if (j < idz) then
              jl = 1
          else 
              jl = idz
          end if
          if (j /= 1 .and. j/= idz .and. abs(zmat(knew,j)) >  zero) then
              call givens(zmat(knew,jl), zmat(knew,j), c, s, r)
              ztemp = zmat(:, j)
              zmat(:, j) = c*zmat(:, j) - s*zmat(:, jl)
              zmat(:, jl) = c*zmat(:, jl) + s*ztemp
              zmat(knew, j) = zero
      !----------------------------------------------------------------!
              !!! In later vesions, we will include the following line.
              !!! For the moment, we exclude it to align with Powell's 
              !!! version.
!-------------!zmat(knew, jl) = r  !--------------------------------! 
      !----------------------------------------------------------------!
          end if
      end do
      
      ! JL plays an important role below. There are two possibilities:
      ! JL = 1 iff IDZ > NPT - N - 1
      ! JL = IDZ iff IDZ <= NPT - N - 1 (but this does not mean JL > 1 
      ! because IDZ may equal 1).

      ! Put the first NPT components of the KNEW-th column of HLAG into 
      ! W, and calculate the parameters of the updating formula.
      tempa = zmat(knew, 1)
      if (idz >=  2) then
          tempa =  - tempa
      end if

      w(1 : npt) = tempa*zmat(:, 1)
      if (jl > 1) then
          tempb = zmat(knew, jl)
          w(1 : npt) = w(1 : npt) + tempb*zmat(:, jl)
      end if

      alpha = w(knew)
      tau = vlag(knew)
      tausq = tau*tau
      denom = alpha*beta + tausq
      vlag(knew) = vlag(knew) - one
      
      ! Complete the updating of ZMAT when there is only one nonzero
      ! element in the KNEW-th row of the new matrix ZMAT, but,  if
      ! IFLAG is set to one, then the first column of ZMAT will be 
      ! exchanged with another one later.
      iflag = 0
      if (jl == 1) then
          temp = sqrt(abs(denom))
          tempb = tempa/temp
          tempa = tau/temp
          zmat(:, 1) = tempa*zmat(:, 1) - tempb*vlag(1 : npt)
          if (idz == 1 .and. temp < zero) then
              ! TEMP < ZERO?!! Powell wrote this but it is STRANGE!!!!!!
              !!! It is possibly a BUG !!!
              idz = 2
          end if
          if (idz >=  2 .and. temp >=  zero) then 
              iflag = 1
          end if
      else
          ! Complete the updating of ZMAT in the alternative case.
          ja = 1
          if (beta >=  zero) then 
              ja = jl
          end if
          jb = jl + 1 - ja
          temp = zmat(knew, jb)/denom
          tempa = temp*beta
          tempb = temp*tau
          temp = zmat(knew, ja)
          scala = one/sqrt(abs(beta)*temp*temp + tausq)
          scalb = scala*sqrt(abs(denom))
          zmat(:, ja) = scala*(tau*zmat(:, ja) - temp*vlag(1 : npt))
          zmat(:, jb) = scalb*(zmat(:, jb) - tempa*w(1 : npt) -         &
     &     tempb*vlag(1 : npt))
          
          if (denom <=  zero) then
              if (beta < zero) then 
                  idz = idz + 1
              end if
              if (beta >=  zero) then 
                  iflag = 1
              end if
          end if
      end if
      
      ! IDZ is reduced in the following case,  and usually the first
      ! column of ZMAT is exchanged with a later one.
      if (iflag == 1) then
          idz = idz - 1
          if (idz > 1) then
              ztemp = zmat(:, 1)
              zmat(:, 1) = zmat(:, idz)
              zmat(:, idz) = ztemp
          end if
      end if
      
      ! Finally,  update the matrix BMAT.
      w(npt + 1 : npt + n) = bmat(knew, :)
      v1 = (alpha*vlag(npt+1 : npt+n) - tau*w(npt+1 : npt+n))/denom
      v2 = (-beta*w(npt+1 : npt+n) - tau*vlag(npt+1 : npt+n))/denom

!      !-------------------POWELL'S IMPLEMENTATION----------------------!
!      do j = 1, n
!          bmat(1 : npt + j, j) = bmat(1 : npt + j, j) +                 &
!     &     v1(j)*vlag(1 : npt + j) + v2(j)*w(1 : npt + j)
!      ! Set the lower triangular part of BMAT(NPT+1:NPT+N,:) by symmetry
!      ! Note that SHIFTBASE sets the upper triangular part by copying
!      ! the lower triangular part, but here it does the opposite. There 
!      ! seems not any particular reason to keep them different. It was
!      ! probably an ad-hoc decision that Powell made when coding. 
!          bmat(npt + j, 1 : j - 1) = bmat(npt + 1 : npt + j - 1, j)
!      end do 
!      !---------------POWELL'S IMPLEMENTATION ENDS---------------------!

      !-----------------MATRIX-VECTOR IMPLEMENTATION-------------------!
      bmat = bmat + outprod(vlag, v1) + outprod(w, v2)
      ! Set the lower triangular part of BMAT(NPT+1:NPT+N,:) by symmetry
      ! Note that SHIFTBASE sets the upper triangular part by copying
      ! the lower triangular part, but here it does the opposite. There 
      ! seems not any particular reason to keep them different. It was
      ! probably an ad-hoc decision that Powell made when coding. 
      ! This part can be spared if we put a pair of parenthsis around 
      ! the two outter products as elaborated below.
      do j = 1, n
          bmat(npt + j, 1 : j - 1) = bmat(npt + 1 : npt + j - 1, j)
      end do 
      !--------------MATRIX-VECTOR IMPLEMENTATION ENDS-----------------!

!---------------------A PROBABLY BETTER IMPLEMENTATION-----------------!
!-----!bmat = bmat + ( outprod(vlag, v1) + outprod(w, v2) ) !----------!
!-----!bmat = bmat + ( outprod(vlag, v1) + outprod(w, v2) ) !----------!
      ! The only difference from the previous matrix-vector 
      ! implementation is the parentheses.
      ! Note that the update is naturally symmetric thanks to the
      ! commutativity of floating point addition. We do not take this
      ! implementation for the moment to produce the same results as
      ! Powell's code, but we should take it in future versions.
!---------------------A PROBABLY BETTER IMPLEMENTATION ENDS------------!
      
      return

      end subroutine update
