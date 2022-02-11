!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of prelim.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun ZHANG (www.zhangzk.net)
! on 12-Feb-2022.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine prelim(calfun, n, npt, x, xl, xu, rhobeg, iprint, maxfu&
     &n, xbase, xpt, fval, gopt, hq, pq, bmat, zmat, ndim, sl, su, nf, k&
     &opt, f, ftarget, xhist, maxxhist, fhist, maxfhist)

      use, non_intrinsic :: consts_mod, only : RP, IK, ZERO, ONE, TWO, H&
     &ALF
      use, non_intrinsic :: evaluate_mod, only : evaluate
      use, non_intrinsic :: history_mod, only : savehist
      use, non_intrinsic :: linalg_mod, only : inprod, matprod, norm
      use, non_intrinsic :: pintrf_mod, only : OBJ

      implicit real(RP) (a - h, o - z)
      implicit integer(IK) (i - n)

      procedure(OBJ) :: calfun
      integer(IK), intent(in) :: maxxhist
      integer(IK), intent(in) :: maxfhist
      integer(IK), intent(in) :: npt
      integer(IK), intent(out) :: nf
      real(RP), intent(in) :: xl(n)
      real(RP), intent(in) :: xu(n)
      real(RP), intent(inout) :: x(n)
      real(RP), intent(out) :: xhist(n, maxxhist)
      real(RP), intent(out) :: fhist(maxfhist)

! Local variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension xbase(n), xpt(npt, n), fval(npt), gopt(n), hq(n * (n + 1&
     &) / 2), pq(npt), bmat(ndim, *), zmat(npt, *), sl(*), su(*)
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonZERO, BOBYQB will change it to its usual value later.
!     NF is maintaned as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!
!     Set some constants.
!
      rhosq = rhobeg * rhobeg
      recip = ONE / rhosq
      np = n + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      almost_infinity = huge(0.0D0) / 2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to ZERO.
!
      do j = 1, n
          xbase(j) = x(j)
          do k = 1, npt
              xpt(k, j) = ZERO
          end do
          do i = 1, ndim
              bmat(i, j) = ZERO
          end do
      end do
      do ih = 1, (n * np) / 2
          hq(ih) = ZERO
      end do
      do k = 1, npt
          pq(k) = ZERO
          do j = 1, npt - np
              zmat(k, j) = ZERO
          end do
      end do
!
!     Begin the initialization procedure. NF becomes ONE more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
      nf = 0
50    nfm = nf
      nfx = nf - n
      nf = nf + 1
      if (nfm <= 2 * n) then
          if (nfm >= 1 .and. nfm <= n) then
              stepa = rhobeg
              if (su(nfm) == ZERO) stepa = -stepa
              xpt(nf, nfm) = stepa
          else if (nfm > n) then
              stepa = xpt(nf - n, nfx)
              stepb = -rhobeg
              if (sl(nfx) == ZERO) stepb = min(TWO * rhobeg, su(nfx))
              if (su(nfx) == ZERO) stepb = max(-TWO * rhobeg, sl(nfx))
              xpt(nf, nfx) = stepb
          end if
      else
          itemp = (nfm - np) / n
          jpt = nfm - itemp * n - n
          ipt = jpt + itemp
          if (ipt > n) then
              itemp = jpt
              jpt = ipt - n
              ipt = itemp
          end if
          xpt(nf, ipt) = xpt(ipt + 1, ipt)
          xpt(nf, jpt) = xpt(jpt + 1, jpt)
      end if
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
      do j = 1, n
          x(j) = min(max(xl(j), xbase(j) + xpt(nf, j)), xu(j))
          if (xpt(nf, j) == sl(j)) x(j) = xl(j)
          if (xpt(nf, j) == su(j)) x(j) = xu(j)
      end do


!-------------------------------------------------------------------!
!call calfun(n, x, f)
      call evaluate(calfun, x, f)
      call savehist(nf, x, xhist, f, fhist)
!-------------------------------------------------------------------!


      fval(nf) = f
      if (nf == 1) then
          fbeg = f
          kopt = 1
      else if (f < fval(kopt)) then
          kopt = nf
      end if
!
!     Set the nonZERO initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
      if (nf <= 2 * n + 1) then
          if (nf >= 2 .and. nf <= n + 1) then
              gopt(nfm) = (f - fbeg) / stepa
              if (npt < nf + n) then
                  bmat(1, nfm) = -ONE / stepa
                  bmat(nf, nfm) = ONE / stepa
                  bmat(npt + nfm, nfm) = -HALF * rhosq
              end if
          else if (nf >= n + 2) then
              ih = (nfx * (nfx + 1)) / 2
              temp = (f - fbeg) / stepb
              diff = stepb - stepa
              hq(ih) = TWO * (temp - gopt(nfx)) / diff
              gopt(nfx) = (gopt(nfx) * stepb - temp * stepa) / diff
              if (stepa * stepb < ZERO) then
                  if (f < fval(nf - n)) then
                      fval(nf) = fval(nf - n)
                      fval(nf - n) = f
                      if (kopt == nf) kopt = nf - n
                      xpt(nf - n, nfx) = stepb
                      xpt(nf, nfx) = stepa
                  end if
              end if
              bmat(1, nfx) = -(stepa + stepb) / (stepa * stepb)
              bmat(nf, nfx) = -HALF / xpt(nf - n, nfx)
              bmat(nf - n, nfx) = -bmat(1, nfx) - bmat(nf, nfx)
              zmat(1, nfx) = sqrt(TWO) / (stepa * stepb)
              zmat(nf, nfx) = sqrt(HALF) / rhosq
              zmat(nf - n, nfx) = -zmat(1, nfx) - zmat(nf, nfx)
          end if
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
      else
          ih = (ipt * (ipt - 1)) / 2 + jpt
          zmat(1, nfx) = recip
          zmat(nf, nfx) = recip
          zmat(ipt + 1, nfx) = -recip
          zmat(jpt + 1, nfx) = -recip
          temp = xpt(nf, ipt) * xpt(nf, jpt)
          hq(ih) = (fbeg - fval(ipt + 1) - fval(jpt + 1) + f) / temp
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     By Tom (on 04-06-2019):
!     If the evaluation returns an NaN or an infinity value, this
!     subroutine is stopped.
      if (f /= f .or. f > almost_infinity) goto 80
!     By Tom (on 04-06-2019):
!     If the target value is reached, stop the algorithm.
      if (f <= ftarget) goto 80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nf < npt .and. nf < maxfun) goto 50
80    nf = min(nf, npt)
      ! nf = npt + 1 at exit of the loop
      end