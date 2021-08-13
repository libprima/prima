!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of trustregion.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 13-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! TRUSTREGION_MOD is a module providing subroutines concerning the
! trust-region iterations.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.
!
! Last Modified: Friday, August 13, 2021 AM12:17:41

      module trustregion_mod

      implicit none
      private
      public :: trsapp, trrad

      contains


      subroutine trsapp(delta, gq, hq, pq, tol, x, xpt, crvmin, s, info)
! TRSAPP finds an approximate solution to the N-dimensional trust
! region subproblem
!
! min <X+S, GQ> + 0.5*<X+S, HESSIAN*(X+S)> s.t. ||S|| <= DELTA
!
! Note that the HESSIAN here is the sum of an explicit part HQ and
! an implicit part (PQ, XPT):
!
! HESSIAN = HQ + sum_K=1^NPT PQ(K)*XPT(:, K)*XPT(:, K)' .
!
! At return, S will be the approximate solution. CRVMIN will be
! set to the least curvature of HESSIAN along the conjugate
! directions that occur, except that it is set to ZERO if S goes
! all the way to the trust region boundary. QRED is the reduction
! of Q achieved by S. INFO is an exit flag:
! INFO = 0: an approximate solution satisfying one of the
! following conditions is found:
! 1. ||G+HS||/||GBEG|| <= TOL,
! 2. ||S|| = DELTA and <S, -(G+HS)> >= (1 - TOL)*||S||*||G+HS||,
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
!
! See Section 5 of the NEWUOA paper.

! Generic modules
      use consts_mod, only : RP, IK, ONE, TWO, HALF, ZERO, PI, DEBUGGING&
     &, SRNLEN
      use debug_mod, only : errstop, verisize
      use infnan_mod, only : is_nan, is_finite
      use lina_mod, only : Ax_plus_y, inprod, matprod

      implicit none

! Inputs
      real(RP), intent(in) :: delta
      real(RP), intent(in) :: gq(:)
      ! GQ(N)
      real(RP), intent(in) :: hq(:, :)
      ! HQ(N, N)
      real(RP), intent(in) :: pq(:)
      ! PQ(NPT)
      real(RP), intent(in) :: tol
      real(RP), intent(in) :: x(:)
      ! X(N)
      real(RP), intent(in) :: xpt(:, :)
      ! XPT(N, NPT)

! Outputs
      integer(IK), intent(out) :: info
      real(RP), intent(out) :: crvmin
      real(RP), intent(out) :: s(:)
      ! S(N)

! Local variables
      integer(IK) :: i
      integer(IK) :: isav
      integer(IK) :: iter
      integer(IK) :: itermax
      integer(IK) :: iu
      integer(IK) :: n
      integer(IK) :: npt
      real(RP) :: alpha
      real(RP) :: angle
      real(RP) :: bstep
      real(RP) :: cf
      real(RP) :: cth
      real(RP) :: d(size(x))
      real(RP) :: dd
      real(RP) :: delsq
      real(RP) :: dg
      real(RP) :: dhd
      real(RP) :: dhs
      real(RP) :: ds
      real(RP) :: g(size(x))
      real(RP) :: gg
      real(RP) :: gg0
      real(RP) :: ggsav
      real(RP) :: hd(size(x))
      real(RP) :: hs(size(x))
      real(RP) :: hx(size(x))
      real(RP) :: qadd
      real(RP) :: qbeg
      real(RP) :: qmin
      real(RP) :: qnew
      real(RP) :: qred
      real(RP) :: qsav
      real(RP) :: quada
      real(RP) :: quadb
      real(RP) :: reduc
      real(RP) :: sg
      real(RP) :: sgk
      real(RP) :: shs
      real(RP) :: sold(size(x))
      real(RP) :: ss
      real(RP) :: sth
      real(RP) :: t
      real(RP) :: unitang
      logical :: twod_search
      character(len=SRNLEN), parameter :: srname = 'TRSAPP'


! Get and verify the sizes.
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUGGING) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(x, n)
          call verisize(gq, n)
          call verisize(hq, n, n)
          call verisize(pq, npt)
          call verisize(s, n)
      end if

      s = ZERO
      crvmin = ZERO
      qred = ZERO
      info = 2
      ! Default exit flag is 2, i.e., itermax is attained

! Prepare for the first line search.
!----------------------------------------------------------------!
!-----!hx = hessmul(hq, pq, xpt, x) !----------------------------!
      hx = Ax_plus_y(hq, x, matprod(xpt, pq * matprod(x, xpt)))
!----------------------------------------------------------------!
      g = gq + hx
      gg = inprod(g, g)
      gg0 = gg
      d = -g
      dd = gg
      ds = ZERO
      ss = ZERO
      hs = ZERO
      delsq = delta * delta
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
! two-dimensional subspace at each iteration being span(S, -G).
      do iter = 1, itermax
! Check whether to exit due to small GG
          if (gg <= (tol**2) * gg0) then
              info = 0
              exit
          end if
! Set BSTEP to the step length such that ||S+BSTEP*D|| = DELTA
          bstep = (delsq - ss) / (ds + sqrt(ds * ds + dd * (delsq - ss))&
     &)
!----------------------------------------------------------------!
!-----!hd = hessmul(hq, pq, xpt, d) !----------------------------!
          hd = Ax_plus_y(hq, d, matprod(xpt, pq * matprod(d, xpt)))
!----------------------------------------------------------------!
          dhd = inprod(d, hd)

! Set the step-length ALPHA and update CRVMIN.
          if (dhd <= ZERO) then
              alpha = bstep
          else
              alpha = min(bstep, gg / dhd)
              if (iter == 1) then
                  crvmin = dhd / dd
              else
                  crvmin = min(crvmin, dhd / dd)
              end if
          end if
! QADD is the reduction of Q due to the new CG step.
          qadd = alpha * (gg - HALF * alpha * dhd)
! QRED is the reduction of Q up to now.
          qred = qred + qadd
! QADD and QRED will be used in the 2D minimization if any.

! Update S, HS, and GG.
          sold = s
          s = s + alpha * d
          ss = inprod(s, s)
          hs = hs + alpha * hd
          ggsav = gg
          ! Gradient norm square before this iteration
          gg = inprod(g + hs, g + hs)
          ! Current gradient norm square
! We may record g+hs for later usage:
! gnew = g + hs
! Note that we should NOT set g = g + hs, because g contains
! the gradient of Q at x.

! Check whether to exit. This should be done after updating HS
! and GG, which will be used for the 2D minimization if any.
          if (alpha >= bstep .or. ss >= delsq) then
! CG path cuts the boundary. Set CRVMIN to 0.
              crvmin = ZERO
! The only possibility that twod_search is true.
              twod_search = .true.
              exit
          end if

! Check whether to exit due to small QADD
          if (qadd <= tol * qred) then
              info = 1
              exit
          end if

! Exit in case of Inf/NaN in S.
          if (.not. is_finite(sum(abs(s)))) then
              s = sold
              info = -1
              exit
          end if

! Prepare for the next CG iteration.
          d = (gg / ggsav) * d - g - hs
          ! CG direction
          dd = inprod(d, d)
          ds = inprod(d, s)
          if (ds <= ZERO) then
! DS is positive in theory.
              info = -1
              exit
          end if
      end do

      if (ss <= ZERO .or. is_nan(ss)) then
! This may occur for ill-conditioned problems due to rounding.
          info = -1
          twod_search = .false.
      end if

      if (twod_search) then
! At least 1 iteration of 2D minimization
          itermax = max(int(1, kind(itermax)), itermax - iter)
      else
          itermax = 0
      end if

! The 2D minimization
      do iter = 1, itermax
          if (gg <= (tol**2) * gg0) then
              info = 0
              exit
          end if
          sg = inprod(s, g)
          shs = inprod(s, hs)
          sgk = sg + shs
          if (sgk / sqrt(gg * delsq) <= tol - ONE) then
              info = 0
              exit
          end if

! Begin the 2D minimization by calculating D and HD and some
! scalar products.
          t = sqrt(delsq * gg - sgk * sgk)
          d = (delsq / t) * (g + hs) - (sgk / t) * s
!----------------------------------------------------------------!
!-----!hd = hessmul(hq, pq, xpt, d) !----------------------------!
          hd = Ax_plus_y(hq, d, matprod(xpt, pq * matprod(d, xpt)))
!----------------------------------------------------------------!
          dg = inprod(d, g)
          dhd = inprod(hd, d)
          dhs = inprod(hd, s)

! Seek the value of the angle that minimizes Q.
          cf = HALF * (shs - dhd)
          qbeg = sg + cf
          qsav = qbeg
          qmin = qbeg
          isav = 0
          iu = 49
          unitang = (TWO * PI) / real(iu + 1, RP)

          do i = 1, iu
              angle = real(i, RP) * unitang
              cth = cos(angle)
              sth = sin(angle)
              qnew = (sg + cf * cth) * cth + (dg + dhs * cth) * sth
              if (qnew < qmin) then
                  qmin = qnew
                  isav = i
                  quada = qsav
              else if (i == isav + 1) then
                  quadb = qnew
              end if
              qsav = qnew
          end do

          if (isav == 0) then
              quada = qnew
          end if
          if (isav == iu) then
              quadb = qbeg
          end if
          if (abs(quada - quadb) > ZERO) then
              quada = quada - qmin
              quadb = quadb - qmin
              angle = HALF * (quada - quadb) / (quada + quadb)
          else
              angle = ZERO
          end if
          angle = unitang * (real(isav, RP) + angle)

! Calculate the new S.
          cth = cos(angle)
          sth = sin(angle)
          reduc = qbeg - (sg + cf * cth) * cth - (dg + dhs * cth) * sth
          sold = s
          s = cth * s + sth * d

! Exit in case of Inf/NaN in S.
          if (.not. is_finite(sum(abs(s)))) then
              s = sold
              info = -1
              exit
          end if

! Calculate HS. Then test for convergence.
          hs = cth * hs + sth * hd
          gg = inprod(g + hs, g + hs)
          qred = qred + reduc
          if (reduc / qred <= tol) then
              info = 1
              exit
          end if
      end do

      end subroutine trsapp


      function trrad(delta, dnorm, eta1, eta2, gamma1, gamma2, ratio)

! Generic module
      use consts_mod, only : RP, HALF

      implicit none

      real(RP) :: trrad
      real(RP), intent(in) :: delta
      ! Current trust region radius
      real(RP), intent(in) :: dnorm
      ! Norm of current trust region step
      real(RP), intent(in) :: eta1
      ! Ratio threshold for contraction
      real(RP), intent(in) :: eta2
      ! Ratio threshold for expansion
      real(RP), intent(in) :: gamma1
      ! Contraction factor
      real(RP), intent(in) :: gamma2
      ! Expansion factor
      real(RP), intent(in) :: ratio
      ! Reduction ratio

      if (ratio <= eta1) then
          trrad = gamma1 * dnorm
      else if (ratio <= eta2) then
          trrad = max(HALF * delta, dnorm)
      else
          trrad = max(HALF * delta, gamma2 * dnorm)
      end if

! For noisy problems, the following may work better.
!if (ratio <= eta1) then
!trrad = gamma1*dnorm
!else if (ratio <= eta2) then  ! Ensure TRRAD >= DELTA
!trrad = delta
!else  ! Ensure TRRAD > DELTA with a constant factor
!trrad = max(delta*(1.0_RP + gamma2)/2.0_RP, gamma2*dnorm)
!end if

      end function trrad


      end module trustregion_mod