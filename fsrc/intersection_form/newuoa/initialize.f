!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of initialize.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 19-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! INITIALIZE_MOD is a module containing subroutines for initializing
! FVAL, XBASE, XPT, GQ, HQ, PQ, IDZ, ZMAT, and BMAT.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.


      module initialize_mod

      implicit none
      private
      public :: initxf, initq, inith


      contains


      subroutine initxf(calfun, iprint, x, rhobeg, ftarget, ij, kopt, nf&
     &, fhist, fval, xbase, xhist, xpt, info)
! INITXF performs the initialization regarding the interpolation
! points and corresponding function values.

! Generic modules
      use consts_mod, only : RP, IK, ZERO, DEBUGGING, SRNLEN
      use debug_mod, only : errstop, verisize
      use info_mod, only : FTARGET_ACHIEVED, NAN_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf
      use output_mod, only : fmssg

! Solver-specific module
      use pintrf_mod, only : FUNEVAL

      implicit none

! Inputs
      procedure(FUNEVAL) :: calfun
      integer(IK), intent(in) :: iprint
      real(RP), intent(in) :: x(:) ! X(N)
      real(RP), intent(in) :: rhobeg
      real(RP), intent(in) :: ftarget

! Outputs
      integer(IK), intent(out) :: info
      integer(IK), intent(out) :: ij(:, :) ! IJ(2, NPT)
      integer(IK), intent(out) :: kopt
      integer(IK), intent(out) :: nf
      real(RP), intent(out) :: fval(:) ! FVAL(NPT)
      real(RP), intent(out) :: fhist(:) ! FHIST(MAXFHIST)
      real(RP), intent(out) :: xbase(:) ! XBASE(N)
      real(RP), intent(out) :: xhist(:, :) ! XHIST(MAXXHIST)
      real(RP), intent(out) :: xpt(:, :) ! XPT(N, NPT)
! Remark on IJ:
! When K > 2*N + 1, all the entries of XPT(:, K) will be zero except
! that the IJ(1, K) and IJ(2, K) entries will be set to RHOBEG or
! -RHOBEG. Consequently, the Hessian of the quadratic model will get a
! possibly nonzero (IJ(1, K), IJ(2, K)) entry.
! Indeed, IJ(:, 1 : 2*N + 1) is never used.

! Intermediate variables
      integer(IK) :: i
      integer(IK) :: itemp
      integer(IK) :: j
      integer(IK) :: k
      integer(IK) :: khist
      integer(IK) :: maxfhist
      integer(IK) :: maxxhist
      integer(IK) :: n
      integer(IK) :: npt
      integer(IK) :: npt_revised
      real(RP) :: f
      real(RP) :: xtemp(size(x))
      logical :: evaluated(size(fval))
      character(len = 6), parameter :: solver= 'NEWUOA'
      character(len = SRNLEN), parameter :: srname = 'INITXF'


! Get and verify the sizes.
      n = int(size(x), kind(n))
      npt = int(size(fval), kind(npt))
      maxfhist = int(size(fhist), kind(maxfhist))
      maxxhist = int(size(xhist, 2), kind(maxxhist))

      if (DEBUGGING) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(X) or SIZE(FVAL) is invalid')
          end if
          if (size(xhist, 1) /= n .and. maxxhist > 0) then
              call errstop(srname, 'XHIST is nonempty but SIZE(XHIST, 1)&
     & /= SIZE(X)')
          end if
          if (maxfhist*maxxhist > 0 .and. maxfhist /= maxxhist) then
              call errstop(srname, 'FHIST and XHIST are both nonempty bu&
     &t SIZE(FHIST) /= SIZE(XHIST, 2)')
          end if
          call verisize(ij, 2_IK, npt)
          call verisize(xbase, n)
          call verisize(xpt, n, npt)
      end if

! At return,
! INFO = 0: initialization finishes normally
! INFO = FTARGET_ACHIEVED: return because f <= ftarget
! INFO = NAN_X: return because x contains NaN
! INFO = NAN_INF_F: return because f is either NaN or +infinity
      info = 0

! We set ij = 1 in case the initialization aborts due to abnormality. If
! we do not do this, ij will be undefined if the initialization aborts.
      ij = 1

! Set XBASE to X.
      xbase = x

! Initialize XPT to ZERO.
      xpt = ZERO

! Begin the initialization procedure. The coordinates of the
! displacement of the next initial interpolation point from XBASE are
! set in XPT(:, .).

! EVALUATED is a boolean array indicating whether the function value of
! the i-th interpolation point has been evaluated. We need it for a
! portable counting of the number of function evaluations, especially
! if the loop is conducted asynchronously. However, the loop here is
! not fully parallelizable if NPT>2N+1, because the definition
! XPT(;, 2N+2:end) depends on FVAL(1:2N+1).
      evaluated = .false.

! NPT_REVISED equals NPT, unless it turns out necessary to return due to
! abnormality (NaN or Inf occurs, or F < FTARGET).
      npt_revised = npt


! Set XPT, FVAL, KOPT, FOPT, and XOPT.

! Set XPT(:, 2 : N + 1).
      do k = 2, min(npt, int(n + 1, kind(npt)))
          xpt(k - 1, k) = rhobeg
      end do
! Set XPT(:, N+2 : NPT)
      do k = int(n + 2, kind(k)), min(npt, int(2*n + 1, kind(npt)))
          xpt(k - n - 1, k) = -rhobeg
      end do

! Set FVAL(1 : NPT) by evaluating F. Totally parallelizable except for
! FMSSG, which outputs messages to the console or files.
      do k = 1, min(npt, int(2*n + 1, kind(npt)))
          xtemp = xpt(:, k) + xbase
          if (any(is_nan(xtemp))) then
              f = sum(xtemp) ! Set F to NaN. It is necessary.
              info = NAN_X
              npt_revised = 0
              exit
          end if
          call calfun(xtemp, f)
          evaluated(k) = .true.
          fval(k) = f

          if (abs(iprint) >= 3) then
              call fmssg(iprint, k, f, xtemp, solver)
          end if
          if (maxfhist >= 1) then
              khist = mod(k - 1_IK, maxfhist) + 1_IK
              fhist(khist) = f
          end if
          if (maxxhist >= 1) then
              khist = mod(k - 1_IK, maxxhist) + 1_IK
              xhist(:, khist) = xtemp
          end if

! Check whether to exit.
          if (f <= ftarget) then
              info = FTARGET_ACHIEVED
              npt_revised = 0
              exit
          end if
          if (is_posinf(f) .or. is_nan(f)) then
              info = NAN_INF_F
              npt_revised = 0
              exit
          end if
      end do

! Set XPT(:, 2*N + 2 : NPT). It depends on FVAL(2 : 2*N + 1).
      do k = int(2*n + 2, kind(k)), npt_revised
! Decide IJ(:, K).  In general, when NPT = (N+1)*(N+2)/2, we can
! set IJ(1 : NPT - (2*N + 1)) to ANY permutation of
! {(I, J) : 1 <= I > J <= N};
! when NPT < (N+1)*(N+2)/2, IJ(1 : NPT - (2*N + 1)) is the
! first NPT - (2*N - 1) elements of such a permutation. Powell took
! the following permutation:
          itemp = int((k - n - 2)/n, kind(itemp))
          j = int(k - (itemp + 1)*n - 1, kind(j))
          i = j + itemp
          if (i > n) then
              itemp = j
              j = i - n
              i = itemp
          end if
          ij(1, k) = i
          ij(2, k) = j

! The following lines set XPT(;, K) to
! XPT(:, I + 1) or XPT(:, I + N + 1)
! +
! XPT(:, J + 1) or XPT(:, J + N + 1),
! depending on the values of FVAL(I + 1), FVAL(I + N + 1),
! FVAL(J + 1), and FVAL(J + N + 1).
!
! This is the only places where the definition
! of XPT(:, 2*N + 2 : NPT) depends on F(2 : 2*N + 1).
! If we set XPT(:, K) to XPT(:, I + 1) + XPT(:, J + 1)
! regardless of FVAL, then the evaluations of FVAL(1 : NPT)
! can be merged, and they are totally parallelizable; this can be
! benificial if the function evaluations are expensive, which is
! likely the case.
          if (fval(i + 1) <= fval(i + n + 1)) then
              xpt(i, k) = xpt(i, i + 1)
          else
              xpt(i, k) = xpt(i, i + n + 1)
          end if
          if (fval(j + 1) <= fval(j + n + 1)) then
              xpt(j, k) = xpt(j, j + 1)
          else
              xpt(j, k) = xpt(j, j + n + 1)
          end if
      end do

! Set FVAL(2*N + 2 : NPT) by evaluating F. Totally parallelizable except
! FMSSG, which outputs messages to the console or files.
      do k = int(2*n + 2, kind(k)), npt_revised
          xtemp = xpt(:, k) + xbase
          if (any(is_nan(xtemp))) then
              f = sum(xtemp) ! Set F to NaN. It is necessary.
              info = NAN_X
              exit
          end if
          call calfun(xtemp, f)
          evaluated(k) = .true.
          fval(k) = f

          if (abs(iprint) >= 3) then
              call fmssg(iprint, k, f, xtemp, solver)
          end if
          if (maxfhist >= 1) then
              khist = mod(k - 1_IK, maxfhist) + 1_IK
              fhist(khist) = f
          end if
          if (maxxhist >= 1) then
              khist = mod(k - 1_IK, maxxhist) + 1_IK
              xhist(:, khist) = xtemp
          end if

! Check whether to exit.
          if (f <= ftarget) then
              info = FTARGET_ACHIEVED
              exit
          end if
          if (is_posinf(f) .or. is_nan(f)) then
              info = NAN_INF_F
              exit
          end if
      end do

! Set NF, KOPT
      nf = int(count(evaluated), kind(nf))
      kopt = int(minloc(fval, dim = 1, mask = evaluated), kind(kopt))

      end subroutine initxf


      subroutine initq(ij, fval, xpt, gq, hq, pq, info)
! INITQ initializes the quadratic model, which is represented by
! (GQ, HQ, PQ) in the following way:
! the gradient of the model at XBASE is GQ;
! the Hessian of the model is
! HQ + sum_{K=1}^NPT PQ(K)*XPT(:, K)*XPT(:, K)'.

! Generic modules
      use consts_mod, only : RP, IK, ZERO, HALF, DEBUGGING, SRNLEN
      use info_mod, only : NAN_MODEL
      use debug_mod, only : errstop, verisize
      use infnan_mod, only : is_nan

      implicit none

! Inputs
      integer(IK), intent(in) :: ij(:, :) ! IJ(2, NPT)
      real(RP), intent(in) :: fval(:) ! FVAL(NPT)
      real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)

! Outputs
      integer(IK), intent(out) :: info
      real(RP), intent(out) :: gq(:) ! GQ(N)
      real(RP), intent(out) :: hq(:, :) ! HQ(N, N)
      real(RP), intent(out) :: pq(:) ! PQ(NPT)

! Intermediate variables
      integer(IK) :: i
      integer(IK) :: j
      integer(IK) :: k
      integer(IK) :: n
      integer(IK) :: npt
      real(RP) :: fbeg
      real(RP) :: fi
      real(RP) :: fj
      real(RP) :: rhobeg
      real(RP) :: xi
      real(RP) :: xj
      character(len = SRNLEN), parameter :: srname = 'INITQ'


! Get and verify the sizes.
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUGGING) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(fval, npt)
          call verisize(ij, 2_IK, npt)
          call verisize(gq, n)
          call verisize(hq, n, n)
          call verisize(pq, npt)
      end if

      gq = ZERO
      hq = ZERO
      pq = ZERO ! We will not update PQ. It is ZERO at return.

      rhobeg = maxval(abs(xpt(:, 2))) ! Read RHOBEG from XPT.
      fbeg = fval(1)

! Set GQ by forward difference.
      gq(1 : n) = (fval(2 : n + 1) - fbeg)/rhobeg
! If possible, revise GQ to central difference.
      k = min(int(npt - n - 1, kind(n)), n)
      gq(1 : k) = HALF*(gq(1 : k) + (fbeg - fval(n+2 : n+1+k))/rhobeg)

! Set the diagonal of HQ by 2nd-order central finite difference.
      do k = 1, min(int(npt - n - 1, kind(n)), n)
          hq(k, k) = ((fval(k + 1) - fbeg)/rhobeg - (fbeg - fval(k + n +&
     & 1))/rhobeg)/rhobeg
      end do
! When NPT > 2*N + 1, set the off-diagonal entries of HQ.
      do k = int(2*n + 2, kind(k)), npt
! I, J, XI, and XJ will be used below.
          i = ij(1, k)
          j = ij(2, k)
          xi = xpt(i, k)
          xj = xpt(j, k)
          if (xi*xpt(i, i + 1) > ZERO) then
              fi = fval(i + 1)
          else
              fi = fval(i + n + 1)
          end if
          if (xj*xpt(j, j + 1) > ZERO) then
              fj = fval(j + 1)
          else
              fj = fval(j + n + 1)
          end if
! With the XI, XJ, FI, and FJ found above, we have
! FVAL(K) = F(XBASE + XI + XJ),
! FI = F(XBASE + XI),
! FJ = F(XBASE + XJ).
! Thus the HQ(I, J) defined below approximates
! frac{partial^2}{partial X_I partial X_J} F(XBASE)
          hq(i, j) = (fbeg - fi - fj + fval(k))/(xi*xj)
          hq(j, i) = hq(i, j)
      end do

      if (any(is_nan(gq)) .or. any(is_nan(hq)) .or. any(is_nan(pq)))then
          info = NAN_MODEL
      else
          info = 0
      end if

      end subroutine initq


      subroutine inith(ij, xpt, idz, bmat, zmat, info)
! INITH initializes BMAT and ZMAT.

! Generic modules
      use consts_mod, only : RP, IK, ZERO, ONE, HALF, DEBUGGING, SRNLEN
      use info_mod, only : NAN_MODEL
      use debug_mod, only : errstop, verisize
      use infnan_mod, only : is_nan

      implicit none

! Inputs
      integer(IK), intent(in) :: ij(:, :) ! IJ(2, NPT)
      real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)

! Outputs
      integer(IK), intent(out) :: idz
      integer(IK), intent(out) :: info
      real(RP), intent(out) :: bmat(:, :) ! BMAT(N, NPT + N)
      real(RP), intent(out) :: zmat(:, :) ! ZMAT(NPT, NPT - N - 1)

! Intermediate variables
      integer(IK) :: i
      integer(IK) :: j
      integer(IK) :: k
      integer(IK) :: n
      integer(IK) :: npt
      real(RP) :: recip
      real(RP) :: reciq
      real(RP) :: rhobeg
      real(RP) :: rhosq
      real(RP) :: xi
      real(RP) :: xj
      character(len = SRNLEN), parameter :: srname = 'INITH'


! Get and verify the sizes.
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUGGING) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(ij, 2_IK, npt)
          call verisize(bmat, n, npt + n)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
      end if

! Set IDZ = 1. It will not be changed in the following.
      idz = 1

! Some values to be used for setting BMAT and ZMAT.
      rhobeg = maxval(abs(xpt(:, 2))) ! Read RHOBEG from XPT.
      rhosq = rhobeg*rhobeg
      recip = ONE/rhosq
      reciq = sqrt(HALF)/rhosq

! Initialize BMAT and ZMAT to ZERO.
      bmat = ZERO
      zmat = ZERO

! Set the nonzero initial elements of BMAT.
! When NPT >= 2*N + 1, this defines BMAT completely;
! When NPT <= 2*N, this defines BMAT(1 : NPT - N - 1, :).
      do k = 1, min(int(npt - n - 1, kind(n)), n)
          bmat(k, k + 1) = HALF/rhobeg
          bmat(k, n + k + 1) = -HALF/rhobeg
      end do

! When NPT <= 2*N, set BMAT(NPT - N : N, :)
      do k = npt - n, n
          bmat(k, 1) = -ONE/rhobeg
          bmat(k, k + 1) = ONE/rhobeg
          bmat(k, npt + k) = -HALF*rhosq
      end do

! Set the nonzero initial elements of ZMAT.
! When NPT <= 2*N + 1, this defines ZMAT completely;
! When NPT > 2*N + 1, this defines ZMAT(:, 1 : N).
      do k = 1, min(int(npt - n - 1, kind(n)), n)
          zmat(1, k) = - reciq - reciq
          zmat(k + 1, k) = reciq
          zmat(k + n + 1, k) = reciq
      end do

! When NPT > 2*N+1, set ZMAT(:, N + 1 : NPT - N - 1).
      do k = int(n + 1, kind(k)), int(npt - n - 1, kind(k))
! I, J, XI, and XJ will be used below.
          i = ij(1, k + n + 1)
          j = ij(2, k + n + 1)
          xi = xpt(i, k + n + 1)
          xj = xpt(j, k + n + 1)
          if (xi < ZERO) then
              i = i + n
          end if
          if (xj < ZERO) then
              j = j + n
          end if
          zmat(1, k) = recip
          zmat(k + n + 1, k) = recip
          zmat(i + 1, k) = -recip
          zmat(j + 1, k) = -recip
      end do

      if (any(is_nan(bmat)) .or. any(is_nan(zmat))) then
          info = NAN_MODEL
      else
          info = 0
      end if

      end subroutine inith


      end module initialize_mod