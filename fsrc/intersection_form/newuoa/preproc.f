!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of preproc.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 07-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! PREPROC_MOD is a module that preprocess the inputs.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.
!
! Last Modified: Saturday, May 22, 2021 PM04:13:55

      module preproc_mod

      implicit none
      private
      public :: preproc


      contains


      subroutine preproc(n, iprint, maxfun, maxhist, npt, eta1, eta2, ft&
     &arget, gamma1, gamma2, rhobeg, rhoend)

      use consts_mod, only : RP, IK, ZERO, ONE, TWO, HALF, TEN, TENTH, E&
     &PS
      use consts_mod, only : RHOBEG_DFT, RHOEND_DFT, FTARGET_DFT, IPRINT&
     &_DFT
      use infnan_mod, only : is_nan, is_inf, is_finite

      implicit none

! Input
      integer(IK), intent(in) :: n

! In-outputs
      integer(IK), intent(inout) :: iprint
      integer(IK), intent(inout) :: maxfun
      integer(IK), intent(inout) :: maxhist
      integer(IK), intent(inout) :: npt
      real(RP), intent(inout) :: eta1
      real(RP), intent(inout) :: eta2
      real(RP), intent(inout) :: ftarget
      real(RP), intent(inout) :: gamma1
      real(RP), intent(inout) :: gamma2
      real(RP), intent(inout) :: rhobeg
      real(RP), intent(inout) :: rhoend

! Local variable
      character(len=6), parameter :: solver = 'NEWUOA'


      if (iprint /= 0 .and. abs(iprint) /= 1 .and. abs(iprint) /= 2 .and&
     &. abs(iprint) /= 3) then
          iprint = IPRINT_DFT
          print '(/1A, I2, 1A)', solver//': invalid IPRINT; it should be&
     & 0, 1, -1, 2, -2, 3, or -3; it is set to ', iprint, '.'
      end if

      if (maxfun < n + 3) then
          maxfun = int(n + 3, kind(maxfun))
          print '(/1A, I8, 1A)', solver//': invalid MAXFUN; it should an&
     & integer at least N + 3 ; it is set to ', maxfun, '.'
      end if

      if (maxhist < 0) then
          maxhist = maxfun
          print '(/1A, I8, 1A)', solver//': invalid MAXHIST; it should b&
     &e a nonnegative integer; it is set to ', maxhist, '.'
      end if
! MAXHIST > MAXFUN is never needed.
      maxhist = min(maxhist, maxfun)

      if (npt < n + 2 .or. npt > min(maxfun - 1, ((n + 2) * (n + 1)) / 2&
     &)) then
          npt = int(min(maxfun - 1, 2 * n + 1), kind(npt))
          print '(/1A, I6, 1A)', solver//': invalid NPT; it should an in&
     &teger in the interval [N+2, (N+1)(N+2)/2], '// 'and it should be l&
     &ess than MAXFUN; it is set to ', npt, '.'
      end if

      if (is_nan(ftarget)) then
          ftarget = FTARGET_DFT
          print '(/1A, 1PD15.6, 1A)', solver//': invalid FTARGET; it sho&
     &uld a real number; it is set to ', ftarget, '.'
      end if

! When the difference between ETA1 and ETA2 is tiny, we force them to equal.
! See the explanation around RHOBEG and RHOEND for the reason.
      if (abs(eta1 - eta2) < 1.0E2_RP * EPS * max(abs(eta1), ONE)) then
          eta2 = eta1
      end if

      if (is_nan(eta1)) then
! In this case, we take the value hard coded in Powell's orginal code
! without any warning. It is useful when intefacing with MATLAB/Python.
          eta1 = TENTH
      else if (eta1 < ZERO .or. eta1 >= ONE) then
! Take ETA1 into account if it has a valid value.
          if (eta2 > ZERO .and. eta2 <= ONE) then
              eta1 = max(EPS, eta2 / 7.0_RP)
          else
              eta1 = TENTH
          end if
          print '(/1A, 1PD15.6, 1A)', solver//': invalid ETA1; it should&
     & be in the interval [0, 1) and not more than ETA2;'// ' it is set &
     &to ', eta1, '.'
      end if

      if (is_nan(eta2)) then
! In this case, we take the value hard coded in Powell's orginal code
! without any warning. It is useful when intefacing with MATLAB/Python.
          eta2 = 0.7_RP
      else if (eta2 < eta1 .or. eta2 > ONE) then
          eta2 = (eta1 + TWO) / 3.0_RP
          print '(/1A, 1PD15.6, 1A)', solver//': invalid ETA2; it should&
     & be in the interval [0, 1] and not less than ETA1;'// ' it is set &
     &to ', eta2, '.'
      end if

      if (is_nan(gamma1)) then
! In this case, we take the value hard coded in Powell's orginal code
! without any warning. It is useful when intefacing with MATLAB/Python.
          gamma1 = HALF
      else if (gamma1 <= ZERO .or. gamma1 >= ONE) then
          gamma1 = HALF
          print '(/1A, 1PD15.6, 1A)', solver//': invalid GAMMA1; it shou&
     &ld in the interval (0, 1); it is set to ', gamma1, '.'
      end if

      if (is_nan(gamma2)) then
! In this case, we take the value hard coded in Powell's orginal code
! without any warning. It is useful when intefacing with MATLAB/Python.
          gamma2 = TWO
      else if (gamma2 < ONE .or. is_inf(gamma2)) then
          gamma2 = TWO
          print '(/1A, 1PD15.6, 1A)', solver//': invalid GAMMA2; it shou&
     &ld a real number not less than 1; it is set to ', gamma2, '.'
      end if

      if (abs(rhobeg - rhoend) < 1.0E2_RP * EPS * max(abs(rhobeg), ONE))&
     & then
! When the data is passed from the interfaces (e.g., MEX) to the Fortran
! code, RHOBEG, and RHOEND may change a bit. It was oberved in a MATLAB
! test that MEX passed 1 to Fortran as 0.99999999999999978. Therefore,
! if we set RHOEND = RHOBEG in the interfaces, then it may happen that
! RHOEND > RHOBEG, which is considered as an invalid input. To avoid this
! situation, we force RHOBEG and RHOEND to equal when the difference is tiny.
          rhoend = rhobeg
      end if

      if (rhobeg <= ZERO .or. is_nan(rhobeg) .or. is_inf(rhobeg)) then
! Take RHOEND into account if it has a valid value.
          if (is_finite(rhoend) .and. rhoend > ZERO) then
              rhobeg = max(TEN * rhoend, RHOBEG_DFT)
          else
              rhobeg = RHOBEG_DFT
          end if
          print '(/1A, 1PD15.6, 1A)', solver//': invalid RHOBEG; it shou&
     &ld be a positive number; it is set to ', rhobeg, '.'
      end if

      if (rhoend <= ZERO .or. rhobeg < rhoend .or. is_nan(rhoend) .or. i&
     &s_inf(rhoend)) then
          rhoend = max(EPS, min(TENTH * rhobeg, RHOEND_DFT))
          print '(/1A, 1PD15.6, 1A)', solver//': invalid RHOEND; it shou&
     &ld be a positive number and RHOEND <= RHOBEG; '// 'it is set to ',&
     & rhoend, '.'
      end if

      end subroutine


      end module preproc_mod