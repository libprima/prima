      subroutine trsapp(n, npt, tol, x, xpt, gq, hq, pq, delta, s,      &
     & crvmin, qred, info)
      ! TRSAPP finds an approximate solution to the N-dimensional trust
      ! region subproblem
      !
      ! min <X+S, GQ> + 0.5*<X+S, HESSIAN*(X+S)> s.t. ||S|| <= DELTA
      !
      ! Note that the HESSIAN here is the sum of an explicit part HQ and
      ! an implicit part (PQ, XPT):
      !
      ! HESSIAN = HQ + sum_K=1^NPT PQ(K)*XPT(K, :)'*XPT(K, :),
      !
      ! where HQ is represented as a vector of its N*(N+1)/2 upper
      ! triangular entries, PQ is an NPT-dimensional vector, and XPT is
      ! an NPT*N matrix.

      ! At return, S will be the approximate solution. CRVMIN will be 
      ! set to the least curvature of HESSIAN along the conjugate 
      ! directions that occur, except that it is set to zero if S goes
      ! all the way to the trust region boundary. QRED is the reduction
      ! of Q achieved by S. INFO is an exit flag:
      ! INFO = 0: an approximate solution satisfying one of the 
      ! following conditions is found: 
      ! 1. ||G||/||GBEG|| <= TOL, 
      ! 2. ||S|| = DELTA and <S, -G> >= (1 - TOL)*||S||*||G||,
      ! where TOL is a tolerance that is set to 1e-2 in NEWUOA.  
      ! INFO = 1: the iteration is reducing Q only slightly;
      ! INFO = 2: the maximal number of iterations is attained;
      ! INFO = -1: too much rounding error to continue
      
      ! The calculation of S begins with the truncated conjugate
      ! gradient method. If the boundary of the trust region is reached,
      ! then further changes to S may be made, each one being in the 2D
      ! space spanned by the current S and the corresponding gradient of
      ! Q. Thus S should provide a substantial reduction to Q within the
      ! trust region.
      
      use consts, only : rp, one, two, half, zero, pi
      use infnan
      use lina
      implicit none
      
      integer, intent(in) :: n, npt
      integer, intent(out) :: info
      real(kind = rp), intent(in) :: tol, delta, x(n), xpt(npt, n),     &
     & gq(n), hq((n*(n + 1))/2), pq(npt)
      real(kind = rp), intent(out) :: s(n), crvmin, qred
      
      integer :: i, isave, iterc, itermax, iu
      real(kind = rp) :: alpha, angle, bstep, cf, cth, dd, delsq, dg,   &
     & dhd, dhs, ds, gg, ggbeg, ggsave, qadd, qbeg, qmin, qnew, qsave,  &
     & reduc, sg, sgk, shs, ss, sth, temp, tempa, tempb, d(n), g(n),    &
     & hd(n), hs(n) 
      logical :: twod_search
      
      s = zero
      crvmin = zero
      qred = zero
      info = 2  ! Default exit flag is 2, i.e., itermax is attained 

      ! Prepare for the first line search.
      call hessmul(n, npt, xpt, hq, pq, x, hd)  ! HD = HESSIAN*X
      g = gq + hd
      gg = dot_product(g, g)
      ggbeg = gg
      d = -g
      dd = gg 
      ds = zero
      ss = zero
      hs = zero
      delsq = delta*delta
      itermax = n

      twod_search = .false. 

      ! The truncated-CG iterations. 
      !
      ! The iteration will be terminated in 4 possible cases:
      ! 1. the maximal number of iterations is attained;
      ! 2. QADD <= TOL*QRED or ||G|| <= TOL*||GBEG||, where QADD is the
      !    reduction of Q due to the latest CG step, QRED is the
      !    reduction of Q since the begnning until the latest CG step,
      !    G is the current gradient, and GBEG is the initial gradient; 
      !    see (5.13) of the NEWUOA paper;
      ! 3. DS <= 0
      ! 4. ||S|| = DELTA, i.e., CG path cuts the trust region boundary.
      ! 
      ! In the 4th case, twod_search will be set to true, meaning that S
      ! will be improved by a sequence of two-dimensional search, the
      ! two-dimensional subspace at each iteration being span(S,-G).
      do iterc = 1, itermax
          ! Check whether to exit due to small GG 
          if (gg <= (tol**2)*ggbeg) then
              info = 0
              exit
          end if
          ! Set BSTEP to the step length such that ||S+BSTEP*D|| = DELTA
          bstep = (delsq-ss)/(ds + sqrt(ds*ds + dd*(delsq-ss)))  
          call hessmul(n, npt, xpt, hq, pq, d, hd)  ! HD = HESSIAN*D
          dhd = dot_product(d, hd)

          ! Set the step-length ALPHA and update CRVMIN and 
          if (dhd <= zero) then
              alpha = bstep
          else
              alpha = min(bstep, gg/dhd)
              if (iterc == 1) then
                  crvmin = dhd/dd
              else
                  crvmin = min(crvmin, dhd/dd)
              end if
          end if
          ! QADD is the reduction of Q due to the new CG step.
          qadd = alpha*(gg - half*alpha*dhd)  
          ! QRED is the reduction of Q up to now.
          qred = qred + qadd
          ! QADD and QRED will be used in the 2D minimization if any.

          ! Update S, HS, and GG.
          s = s + alpha*d 
          ss = dot_product(s, s)
          hs = hs + alpha*hd
          ggsave = gg  ! Gradient norm square before this iteration 
          gg = dot_product(g+hs, g+hs)  ! Current gradient norm square
             
          ! Check whether to exit. This should be done after updating HS
          ! and GG, which will be used for the 2D minimization if any. 
          if (alpha >= bstep .or. ss >= delsq) then
              ! CG path cuts the boundary. Set CRVMIN to 0.
              crvmin = zero
              ! The only possibility that twod_search is true.
              twod_search = .true.
              exit
          end if
           
          ! Check whether to exit due to small QADD
          if (qadd <= tol*qred) then
              info = 1 
              exit
          end if

          ! Prepare for the next CG iteration.
          d = (gg/ggsave)*d - g - hs  ! CG direction
          dd = dot_product(d, d)
          ds = dot_product(d, s)
          if (ds <= zero) then 
              ! DS is positive in theory.
              info = -1 
              exit 
          end if
      end do 

      if (ss <= 0 .or. is_nan(ss)) then 
          ! This may occur for ill-conditioned problems due to rounding.
          info = -1
          twod_search = .false.
      end if

      if (twod_search) then
          ! At least 1 iteration of 2D minimization
          itermax = max(1, itermax - iterc)  
      else
          itermax = 0
      end if
      
      ! The 2D minimization 
      do iterc = 1, itermax
          if (gg <= (tol**2)*ggbeg) then 
              info = 0
              exit
          end if
          sg = dot_product(s, g)
          shs = dot_product(s, hs)
          sgk = sg + shs
          if (sgk/sqrt(gg*delsq) <= tol - one) then 
              info = 0
              exit
          end if
          
          ! Begin the 2D minimization by calculating D and HD and some
          ! scalar products.
          temp = sqrt(delsq*gg - sgk*sgk)
          d = (delsq/temp)*(g + hs) - (sgk/temp)*s
          call hessmul(n, npt, xpt, hq, pq, d, hd)  ! HD = HESSIAN*D
          dg = dot_product(d, g)
          dhd = dot_product(hd, d)
          dhs = dot_product(hd, s)

          ! Seek the value of the angle that minimizes Q.
          cf = half*(shs - dhd)
          qbeg = sg + cf
          qsave = qbeg
          qmin = qbeg
          isave = 0
          iu = 49
          temp = two*pi/real(iu + 1, rp)
          do i = 1, iu
              angle = real(i, rp)*temp
              cth = cos(angle)
              sth = sin(angle)
              qnew = (sg + cf*cth)*cth + (dg + dhs*cth)*sth
              if (qnew < qmin) then
                  qmin = qnew
                  isave = i
                  tempa = qsave
              else if (i == isave + 1) then
                  tempb = qnew
              end if
              qsave = qnew
          end do
          if (isave == 0) then
              tempa = qnew
          end if
          if (isave == iu) then 
              tempb = qbeg
          end if
          tempa = tempa - qmin
          tempb = tempb - qmin
          if (abs(tempa - tempb) > zero) then
              angle = half*(tempa - tempb)/(tempa + tempb)
          else
              angle = zero
          end if
          angle = temp*(real(isave, rp) + angle)

          ! Calculate the new S and HS. Then test for convergence.
          cth = cos(angle)
          sth = sin(angle)
          reduc = qbeg - (sg + cf*cth)*cth - (dg + dhs*cth)*sth
          s = cth*s + sth*d
          hs = cth*hs + sth*hd
          gg = dot_product(g+hs, g+hs)
          qred = qred + reduc
          if (reduc/qred <= tol) then 
              info = 1
              exit
          end if
      end do

      return

      end subroutine trsapp
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hessmul(n, npt, xpt, hq, pq, d, hd)
      ! HESSMUL calculates HD = HESSIAN*D for the NEWUOA-form HESSIAN,
      ! which is the sum of an explicit part HQ and an implicit 
      ! part (PQ, XPT). Specifically, 
      !
      ! HESSIAN = HQ + sum_K=1^NPT PQ(K)*XPT(K, :)'*XPT(K, :),
      !
      ! where HQ is represented as a vector of its N*(N+1)/2 upper
      ! triangular entries, PQ is an NPT-dimensional vector, and XPT is
      ! an NPT*N matrix.
      
      use consts, only : rp, zero

      integer, intent(in) :: n, npt
      real(kind = rp), intent(in) :: xpt(npt, n), hq((n*(n + 1))/2),    &
     & pq(npt), d(n)
      
      real(kind = rp), intent(out) :: hd(n)

      integer :: i, ih, j, k

      hd = zero
      do k = 1, npt
          hd = hd + dot_product(xpt(k, :), d)*pq(k)*xpt(k, :)
      end do

      ih = 0
      do j = 1, n
         do i = 1, j
             ih = ih + 1
             if (i < j) then 
                 hd(j) = hd(j) + hq(ih)*d(i)
             end if
             hd(i) = hd(i) + hq(ih)*d(j)
         end do
      end do

      end subroutine
