!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of vlagbeta.f90.
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


      module vlagbeta_mod
!--------------------------------------------------------------------------------------------------!
! This module contains a subroutine that calculates VLAG and BETA for a given step D. Both VLAG and
! BETA are critical for the updating procedure of H, which is detailed formula (4.11) of the NEWUOA
! paper. See (4.12) for the definition of BETA, and VLAG is indeed Hw.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the NEWUOA paper.
!
! Started: July 2020
!
! Last Modified: Friday, November 05, 2021 PM02:02:49
!--------------------------------------------------------------------------------------------------!

      implicit none
      private
      public :: calvlag, calbeta


      contains


      function calvlag(idz, kopt, bmat, d, xpt, zmat) result(vlag)
!--------------------------------------------------------------------------------------------------!
! This function calculates VLAG = Hw for a given step D. See (4.11) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
      use, non_intrinsic :: consts_mod, only : RP, IK, ONE, HALF, DEBUGG&
     &ING
      use, non_intrinsic :: debug_mod, only : assert
      use, non_intrinsic :: infnan_mod, only : is_finite
      use, non_intrinsic :: linalg_mod, only : matprod, omega_mul, issym&
     &metric

      implicit none

! Inputs
      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: kopt
      real(RP), intent(in) :: bmat(:, :)
      ! BMAT(N, NPT + N)
      real(RP), intent(in) :: d(:)
      ! D(N)
      real(RP), intent(in) :: xpt(:, :)
      ! XPT(N, NPT)
      real(RP), intent(in) :: zmat(:, :)
      ! ZMAT(NPT, NPT - N - 1)

! Outputs
      real(RP) :: vlag(size(bmat, 2))
      ! VLAG(NPT + N)

! Local variables
      character(len=*), parameter :: srname = 'CALVLAG'
      integer(IK) :: n
      integer(IK) :: npt
      real(RP) :: wcheck(size(zmat, 1))
      real(RP) :: xopt(size(xpt, 1))

! Sizes
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

! Preconditions
      if (DEBUGGING) then
          call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2',&
     & srname)
          call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ&
     & <= SIZE(ZMAT, 2) + 1', srname)
          call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', s&
     &rname)
          call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n,&
     & 'SIZE(BMAT)==[N, NPT+N]', srname)
          call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NP&
     &T+1:NPT+N) is symmetric', srname)
          call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - &
     &n - 1, 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
          call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == &
     &N, D is finite', srname)
          call assert(all(is_finite(xpt)), 'XPT is finite', srname)
      end if

!====================!
! Calculation starts !
!====================!

      xopt = xpt(:, kopt)
      ! Read XOPT.

! WCHECK contains the first NPT entries of (w-v) for the vectors w and v defined in (4.10) and
! (4.24) of the NEWUOA paper, and also hat{w} in eq(6.5) of
! M. J. D. Powell, Least Frobenius norm updating of quadratic models that satisfy interpolation
! conditions. Math. Program., 100:183--215, 2004
      wcheck = matprod(d, xpt)
      wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))

      vlag(1:npt) = matprod(d, bmat(:, 1:npt))
      vlag(1:npt) = vlag(1:npt) + omega_mul(idz, zmat, wcheck)
      vlag(kopt) = vlag(kopt) + ONE
      ! The calculation of VLAG(1:NPT) finishes.
      vlag(npt + 1:npt + n) = matprod(bmat, [wcheck, d])
      ! The calculation of VLAG finishes.

!====================!
!  Calculation ends  !
!====================!

! Postconditions
      if (DEBUGGING) then
          call assert(size(vlag) == npt + n, 'SIZE(VLAG) == NPT + N', sr&
     &name)
      end if

      end function calvlag


      function calbeta(idz, kopt, bmat, d, xpt, zmat) result(beta)
!--------------------------------------------------------------------------------------------------!
! This function calculates BETA for a given step D. See (4.12) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!

! Generic modules
      use, non_intrinsic :: consts_mod, only : RP, IK, TWO, HALF, DEBUGG&
     &ING
      use, non_intrinsic :: debug_mod, only : assert
      use, non_intrinsic :: infnan_mod, only : is_finite
      use, non_intrinsic :: linalg_mod, only : inprod, matprod, omega_in&
     &prod, issymmetric

      implicit none

! Inputs
      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: kopt
      real(RP), intent(in) :: bmat(:, :)
      ! BMAT(N, NPT + N)
      real(RP), intent(in) :: d(:)
      ! D(N)
      real(RP), intent(in) :: xpt(:, :)
      ! XPT(N, NPT)
      real(RP), intent(in) :: zmat(:, :)
      ! ZMAT(NPT, NPT - N - 1)

! Outputs
      real(RP) :: beta

! Local variables
      character(len=*), parameter :: srname = 'CALBETA'
      integer(IK) :: n
      integer(IK) :: npt
      real(RP) :: bw(size(bmat, 1))
      real(RP) :: bd(size(bmat, 1))
      real(RP) :: bsum
      real(RP) :: dsq
      real(RP) :: dx
      real(RP) :: wcheck(size(zmat, 1))
      real(RP) :: xopt(size(xpt, 1))
      real(RP) :: xoptsq

! Sizes
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

! Preconditions
      if (DEBUGGING) then
          call assert(n >= 1 .and. npt >= n + 2, 'N >= 1, NPT >= N + 2',&
     & srname)
          call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ&
     & <= SIZE(ZMAT, 2) + 1', srname)
          call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', s&
     &rname)
          call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n,&
     & 'SIZE(BMAT)==[N, NPT+N]', srname)
          call assert(issymmetric(bmat(:, npt + 1:npt + n)), 'BMAT(:, NP&
     &T+1:NPT+N) is symmetric', srname)
          call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - &
     &n - 1, 'SIZE(ZMAT) == [NPT, NPT - N - 1]', srname)
          call assert(size(d) == n .and. all(is_finite(d)), 'SIZE(D) == &
     &N, D is finite', srname)
          call assert(all(is_finite(xpt)), 'XPT is finite', srname)
      end if

!====================!
! Calculation starts !
!====================!

      xopt = xpt(:, kopt)
      ! Read XOPT.

      dx = inprod(d, xopt)
      dsq = inprod(d, d)
      xoptsq = inprod(xopt, xopt)

! WCHECK contains the first NPT entries of (w-v) for the vectors w and v defined in (4.10) and
! (4.24) of the NEWUOA paper, and also hat{w} in eq(6.5) of
! M. J. D. Powell, Least Frobenius norm updating of quadratic models that satisfy interpolation
! conditions. Math. Program., 100:183--215, 2004
      wcheck = matprod(d, xpt)
      wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))

      bw = matprod(bmat(:, 1:npt), wcheck)
      bd = matprod(bmat(:, npt + 1:npt + n), d)
      bsum = sum(bd * d + bw * d + bw * d)

      beta = dx**2 + dsq * (xoptsq + TWO * dx + HALF * dsq) - omega_inpr&
     &od(idz, zmat, wcheck, wcheck) - bsum

!====================!
!  Calculation ends  !
!====================!

      end function calbeta

      end module vlagbeta_mod