      module update_mod

      implicit none
      private
      public :: updateh, updateq, qalt


      contains

      subroutine updateh(bmat, zmat, idz, vlag, beta, knew)
      ! UPDATE updates arrays BMAT and ZMAT together with IDZ, in order
      ! to shift the interpolation point that has index KNEW. On entry, 
      ! VLAG contains the components of the vector THETA*WCHECK + e_b 
      ! of the updating formula (6.11) in the NEWUOA paper, and BETA
      ! holds the value of the parameter that has this name. 
      
      ! Although VLAG will be modified below, its value will NOT be used 
      ! by other parts of NEWUOA after returning from UPDATE. Its value 
      ! will be overwritten when trying the alternative model or by
      ! VLAGBETA.

      use consts_mod, only : RP, IK, ONE, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none

      integer(IK), intent(in) :: knew
      integer(IK), intent(inout) :: idz
      real(RP), intent(in) :: beta
      real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
      real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
      real(RP), intent(inout) :: vlag(:)     ! VLAG(NPT + N)

      integer(IK) :: iflag, j, ja, jb, jl, n, npt
      real(RP) :: alpha, denom, scala, scalb, tau, tausq, temp,         &
     & tempa, tempb, ztemp(size(zmat, 1)), w(size(vlag)),               &
     & v1(size(bmat, 1)), v2(size(bmat, 1))
      character(len = SRNLEN), parameter :: srname = 'UPDATEH'

      
      ! Get and verify the sizes.
      n = int(size(bmat, 1), kind(n))
      npt = int(size(bmat, 2), kind(npt)) - n

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(BMAT) is invalid')
          end if
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(vlag, npt + n)
      end if
      
          
      ! Apply the rotations that put zeros in the KNEW-th row of ZMAT.
      ! A Givens rotation will be multiplied to ZMAT from the left so
      ! ZMAT(KNEW, JL) becomes SQRT(ZMAT(KNEW, JL)^2+ZMAT(KNEW,J)) and
      ! ZMAT(KNEW, J) becomes 0. 
      jl = 1  ! For J = 2, ..., IDZ - 1, set JL = 1.
      do j = 2, int(idz - 1, kind(j))
          call grota(zmat, jl, j, knew)
      end do

      if (idz <= npt - n - 1) then
          jl = idz  ! For J = IDZ + 1, ..., NPT - N - 1, set JL = IDZ.
      end if
      do j = int(idz + 1, kind(j)), int(npt - n - 1, kind(j))
          call grota(zmat, jl, j, knew)
      end do
      
      ! JL plays an important role below. There are two possibilities:
      ! JL = 1 < IDZ iff IDZ = 1 
      ! JL = IDZ > 1 iff 2 <= IDZ <= NPT - N - 1

      ! Put the first NPT components of the KNEW-th column of HLAG into 
      ! W, and calculate the parameters of the updating formula.
      tempa = zmat(knew, 1)
      if (idz >=  2) then
          tempa = -tempa
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
      vlag(knew) = vlag(knew) - ONE
      
      ! Complete the updating of ZMAT when there is only one nonzero
      ! element in the KNEW-th row of the new matrix ZMAT, but, if
      ! IFLAG is set to one, then the first column of ZMAT will be 
      ! exchanged with another one later.
      iflag = 0
      if (jl == 1) then
          ! There is only one nonzero in ZMAT(KNEW, :) after the
          ! rotation. This is the normal case, because IDZ is 1 in
          ! precise arithmetic. 
          temp = sqrt(abs(denom))
          tempb = tempa/temp
          tempa = tau/temp
          zmat(:, 1) = tempa*zmat(:, 1) - tempb*vlag(1 : npt)
      !----------------------------------------------------------------!
          if (idz == 1 .and. temp < ZERO) then
!---------!if (idz == 1 .and. denom < ZERO) then !------------------!
      !----------------------------------------------------------------!
              ! TEMP < ZERO?!! Powell wrote this but it is STRANGE!!!!!!
              !!! It is probably a BUG !!!
              ! According to (4.18) of the NEWUOA paper, the 
              ! "TEMP < ZERO" here and "TEMP >= ZERO" below should be
              ! revised by replacing "TEMP" with DENOM, which is denoted
              ! by sigma in the paper. See also the corresponding part
              ! of the LINCOA code (which has also some strangeness).
              ! It seems that the BOBYQA code does not have this part
              ! --- it does not have IDZ at all (why?).
              idz = 2
          end if
      !----------------------------------------------------------------!
          if (idz >= 2 .and. temp >= ZERO) then 
!---------!if (idz >= 2 .and. denom >= ZERO) then !--------------------!
      !----------------------------------------------------------------!
              ! JL = 1 and IDZ >= 2??? Seems not possible either!!!
              iflag = 1
          end if
      else
          ! Complete the updating of ZMAT in the alternative case.
          ! There are two nonzeros in ZMAT(KNEW, :) after the rotation.
          ja = 1
          if (beta >=  ZERO) then 
              ja = jl
          end if
          jb = int(jl + 1 - ja, kind(jb))
          temp = zmat(knew, jb)/denom
          tempa = temp*beta
          tempb = temp*tau
          temp = zmat(knew, ja)
          scala = ONE/sqrt(abs(beta)*temp*temp + tausq)
          scalb = scala*sqrt(abs(denom))
          zmat(:, ja) = scala*(tau*zmat(:, ja) - temp*vlag(1 : npt))
          zmat(:, jb) = scalb*(zmat(:, jb) - tempa*w(1 : npt) -         &
     &     tempb*vlag(1 : npt))
          
          if (denom <=  ZERO) then
              if (beta < ZERO) then 
                  idz = int(idz + 1, kind(idz))  
                  ! Is it possible to get IDZ>NPT-N-1?
              end if
              if (beta >=  ZERO) then 
                  iflag = 1
              end if
          end if
      end if
      
      ! IDZ is reduced in the following case,  and usually the first
      ! column of ZMAT is exchanged with a later one.
      if (iflag == 1) then
          idz = int(idz - 1, kind(idz))
          if (idz > 1) then
              ztemp = zmat(:, 1)
              zmat(:, 1) = zmat(:, idz)
              zmat(:, idz) = ztemp
          end if
      end if
      
      ! Finally,  update the matrix BMAT.
      w(npt + 1 : npt + n) = bmat(:, knew)
      v1 = (alpha*vlag(npt+1 : npt+n) - tau*w(npt+1 : npt+n))/denom
      v2 = (-beta*w(npt+1 : npt+n) - tau*vlag(npt+1 : npt+n))/denom

      call r2update(bmat, ONE, v1, vlag, ONE, v2, w)
      ! In floating-point arithmetic, the update above does not guarante
      ! BMAT(:, NPT+1 : NPT+N) to be symmetric. Symmetrization needed.
      call symmetrize(bmat(:, npt + 1 : npt + n))

      return

      end subroutine updateh


      subroutine updateq(idz, knew, fqdiff, xptknew, bmatknew, zmat, gq,&
     & hq, pq)

      use warnerror_mod, only : errstop
      use consts_mod, only : RP, IK, ZERO, DEBUG_MODE, SRNLEN
      use lina_mod
      implicit none

      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: knew
      real(RP), intent(in) :: fqdiff
      real(RP), intent(in) :: xptknew(:)  ! XPTKNEW(N)
      real(RP), intent(in) :: bmatknew(:) ! BMATKNEW(N)
      real(RP), intent(in) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)
      real(RP), intent(inout) :: gq(:)   ! GQ(N)
      real(RP), intent(inout) :: hq(:, :)! HQ(N, N)
      real(RP), intent(inout) :: pq(:)   ! PQ(NPT) 

      integer(IK) :: n, npt
      real(RP) :: fqdz(size(zmat, 2))
      character(len = SRNLEN), parameter :: srname = 'UPDATEQ'

      
      ! Get and verify the sizes.
      n = int(size(gq), kind(n))
      npt = int(size(pq), kind(npt))

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(GQ) or SIZE(PQ) is invalid')
          end if
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(xptknew, n)
          call verisize(bmatknew, n)
          call verisize(hq, n, n)
      end if

      !----------------------------------------------------------------!
      ! Implement R1UPDATE properly so that it ensures HQ is symmetric.
      call r1update(hq, pq(knew), xptknew)
      !----------------------------------------------------------------!

      ! Update the implicit part of second derivatives.
      fqdz = fqdiff*zmat(knew, :)
      fqdz(1 : idz - 1) = -fqdz(1 : idz - 1)
      pq(knew) = ZERO
      !----------------------------------------------------------------!
!----!pq = pq + matmul(zmat, fqdz) !-----------------------------------!
      pq = Ax_plus_y(zmat, fqdz, pq)
      !----------------------------------------------------------------!

      ! Update the gradient.
      gq = gq + fqdiff*bmatknew

      return 

      end subroutine updateq


      subroutine qalt(gq, hq, pq, fval, smat, zmat, kopt, idz)
      ! QALT calculates the alternative model, namely the model that
      ! minimizes the F-norm of the Hessian subject to the interpolation
      ! conditions. 
      ! Note that SMAT = BMAT(:, 1:NPT)

      use consts_mod, only : RP, IK, ZERO, DEBUG_MODE, SRNLEN
      use warnerror_mod, only : errstop
      use lina_mod
      implicit none

      integer(IK), intent(in) :: kopt
      integer(IK), intent(in) :: idz
      real(RP), intent(in) :: fval(:)       ! FVAL(NPT)
      real(RP), intent(in) :: smat(:, :)    ! SMAT(N, NPT)
      real(RP), intent(in) :: zmat(:, :)    ! ZMAT(NPT, NPT-N-!)
      real(RP), intent(out) :: gq(:)        ! GQ(N)
      real(RP), intent(out) :: hq(:, :)     ! HQ(N, N)
      real(RP), intent(out) :: pq(:)        ! PQ(NPT)

      real(RP) :: vlag(size(pq)), vz(size(zmat, 2))
      integer(IK) :: n, npt
      character(len = SRNLEN), parameter :: srname = 'QALT'


      ! Get and verify the sizes.
      n = int(size(gq), kind(n))
      npt = int(size(pq), kind(npt))

      if (DEBUG_MODE) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(GQ) or SIZE(PQ) is invalid')
          end if
          call verisize(fval, npt)
          call verisize(smat, n, npt)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(hq, n, n)
      end if

      vlag = fval - fval(kopt)
      gq = matmul(smat, vlag)
      hq = ZERO
      vz = matmul(vlag, zmat)
      vz(1 : idz - 1) = - vz(1 : idz - 1)
      pq = matmul(zmat, vz)

      return

      end subroutine qalt

      end module update_mod
