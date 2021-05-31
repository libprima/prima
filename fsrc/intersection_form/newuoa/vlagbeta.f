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
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 31-May-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! VLAGBETA_MOD is a module providing a suroutine that calculates VLAG
! and BETA for a given step D.
! Both VLAG and BETA are critical for the updating procedure of H, which
! is detailed formula (4.11) of the NEWUOA paper. See (4.12) for the
! definition of BETA, and VLAG is indeed Hw.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.
!
! Last Modified: Saturday, May 22, 2021 PM04:18:51

      module vlagbeta_mod

      implicit none
      private
      public :: vlagbeta


      contains


      subroutine vlagbeta(idz, kopt, bmat, d, xopt, xpt, zmat, beta, vla&
     &g)
! VLAGBETA is calculates VLAG = Hw and BETA for a given step D.
! See (4.11)--(4.12) of the NEWUOA paper.

! Generic modules
      use consts_mod, only : RP, IK, ONE, HALF, DEBUGGING, SRNLEN
      use debug_mod, only : errstop, verisize
      use lina_mod, only : Ax_plus_y, xA_plus_y, xpy_dot_z, inprod, matp&
     &rod

      implicit none

! Inputs
      integer(IK), intent(in) :: idz
      integer(IK), intent(in) :: kopt
      real(RP), intent(in) :: bmat(:, :) ! BMAT(N, NPT + N)
      real(RP), intent(in) :: d(:) ! D(N)
      real(RP), intent(in) :: xopt(:) ! XOPT(N)
      real(RP), intent(in) :: xpt(:, :) ! XPT(N, NPT)
      real(RP), intent(in) :: zmat(:, :) ! ZMAT(NPT, NPT - N - 1)

! Outputs
      real(RP), intent(out) :: beta
      real(RP), intent(out) :: vlag(:) ! VLAG(NPT + N)

! Intermediate variables
      integer(IK) :: n
      integer(IK) :: npt
      real(RP) :: bw(size(bmat, 1))
      real(RP) :: bwvd
      real(RP) :: dsq
      real(RP) :: dx
      real(RP) :: wcheck(size(zmat, 1))
      real(RP) :: wz(size(zmat, 2))
      real(RP) :: wzsave(size(wz))
      real(RP) :: xoptsq
      character(len=SRNLEN), parameter :: srname = 'VLAGBETA'


! Get and verify the sizes
      n = int(size(xpt, 1), kind(n))
      npt = int(size(xpt, 2), kind(npt))

      if (DEBUGGING) then
          if (n == 0 .or. npt < n + 2) then
              call errstop(srname, 'SIZE(XPT) is invalid')
          end if
          call verisize(bmat, n, npt + n)
          call verisize(zmat, npt, int(npt - n - 1, kind(n)))
          call verisize(xopt, n)
          call verisize(d, n)
          call verisize(vlag, n + npt)
      end if


!----------------------------------------------------------------------!
! This is the one of the two places where WCHECK is calculated,
! the other one being BIGDEN (now removed).
! WCHECK contains the first NPT entries of (w-v) for the vectors
! w and v defined in eq(4.10) and eq(4.24) of the NEWUOA paper,
! and also hat{w} in eq(6.5) of
!
! M. J. D. Powell, Least Frobenius norm updating of quadratic
! models that satisfy interpolation conditions. Math. Program.,
! 100:183--215, 2004
      wcheck = matprod(d, xpt)
      wcheck = wcheck * (HALF * wcheck + matprod(xopt, xpt))
!----------------------------------------------------------------------!

      vlag(1:npt) = matprod(d, bmat(:, 1:npt))

      wz = matprod(wcheck, zmat)
      wzsave = wz
      wz(1:idz - 1) = -wz(1:idz - 1)
      beta = -inprod(wzsave, wz)
!----------------------------------------------------------------------!
!-----!vlag(1 : npt) = vlag(1 : npt) + matprod(zmat, wz) !-------------!
      vlag(1:npt) = Ax_plus_y(zmat, wz, vlag(1:npt))
!----------------------------------------------------------------------!

      bw = matprod(bmat(:, 1:npt), wcheck)
!----------------------------------------------------------------------!
!vlag(npt + 1 : npt + n) = bw + matprod(d, bmat(:, npt + 1 : npt + n)) !
      vlag(npt + 1:npt + n) = xA_plus_y(bmat(:, npt + 1:npt + n), d, bw)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!-----!bwvd = inprod(bw + vlag(npt + 1 : npt + n), d) !----------------!
      bwvd = xpy_dot_z(bw, vlag(npt + 1:npt + n), d)
!----------------------------------------------------------------------!

      dx = inprod(d, xopt)

      dsq = inprod(d, d)
      xoptsq = inprod(xopt, xopt)

      beta = dx * dx + dsq * (xoptsq + dx + dx + HALF * dsq) + beta - bw&
     &vd
      vlag(kopt) = vlag(kopt) + ONE

      end subroutine vlagbeta


      end module vlagbeta_mod