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
! on 11-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
     &_DFT, MAXIMAL_HIST
      use infnan_mod, only : is_nan, is_inf

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

! Intermediate variables
      character(len = 6), parameter :: solver = 'NEWUOA'


      if (iprint /= 0 .and. iprint /= 1 .and. iprint /= 2 .and. iprint /&
     &= 3 .and. iprint /= 4) then
          iprint = IPRINT_DFT
          print '(/1X, 1A, I1, 1A)', solver // ': invalid IPRINT; it sho&
     &uld be 0, 1, 2, 3, or 4; it is set to ', iprint, '.'
      end if

      if (maxfun < n + 3) then
          maxfun = int(n + 3, kind(n))
          print '(/1X, 1A, I6, 1A)', solver // ': invalid MAXFUN; it sho&
     &uld an integer at least N + 3 ; it is set to ', maxfun, '.'
      end if

      if (maxhist < 0 .or. maxhist > MAXIMAL_HIST) then
          maxhist = max(0_IK, min(MAXIMAL_HIST, maxfun))
          print '(/1X, 1A, I8, 1A, I8, 1A)', solver // ': invalid MAXHIS&
     &T; it should be a nonnegative integer not more thant ', MAXIMAL_HI&
     &ST, ', it is set to ', maxhist, '.'
      end if

      if (npt < n + 2 .or. npt > min(maxfun - 1, ((n + 2)*(n + 1))/2)) t&
     &hen
          npt = int(min(maxfun - 1, 2*n + 1), kind(npt))
          print '(/1X, 1A, I6, 1A)', solver // ': invalid NPT; it should&
     & an integer in the interval [N+2, (N+1)(N+2)/2], ' // 'and it shou&
     &ld be less than MAXFUN; it is set to ', npt, '.'
      end if

      if (eta1 < ZERO .or. eta1 > ONE .or. is_nan(eta1)) then
          eta1 = MAX(EPS, eta2/7.0_RP)
          print '(/1X, 1A, 1PD15.6, 1A)', solver // ': invalid ETA1; it &
     &should be in the interval (0, 1); it is set to ', eta1, '.'
      end if

      if (eta2 < eta1 .or. eta2 > ONE .or. is_nan(eta2)) then
          eta2 = (eta1 + TWO)/3.0_RP
          print '(/1X, 1A, 1PD15.6, 1A)', solver // ': invalid ETA2; it &
     &should in the intercal [ETA1, 1); it is set to ', eta2, '.'
      end if

      if (is_nan(ftarget)) then
          ftarget = FTARGET_DFT
          print '(/1X, 1A, 1PD15.6, 1A)', solver // ': invalid FTARGET; &
     &it should a real number; it is set to ', ftarget, '.'
      end if

      if (gamma1 <= ZERO .or. gamma1 >= ONE .or. is_nan(gamma1)) then
          gamma1 = HALF
          print '(/1X, 1A, 1PD15.6, 1A)', solver // ': invalid GAMMA1; i&
     &t should in the intercal (0, 1); it is set to ', gamma1, '.'
      end if

      if (gamma2 < ONE .or. is_nan(gamma2) .or. is_inf(gamma2)) then
          gamma2 = TWO
          print '(/1X, 1A, 1PD15.6, 1A)', solver // ': invalid GAMMA2; i&
     &t should a real number not less than 1; it is set to ', gamma2, '.&
     &'
      end if

      if ((rhobeg - rhoend) < 1.0e2_RP*EPS*max(abs(rhobeg), ONE))then
! When the data is passed from the interfaces (e.g., MEX) to the Fortran
! code, RHOBEG, and RHOEND may change a bit. It was oberved in a MATLAB
! test that MEX passed 1 to Fortran as 0.99999999999999978. Therefore,
! if we set RHOEND = RHOBEG in the interfaces, then it may happen that
! RHOEND > RHOBEG, which is considered as an invalid input. To avoid this
! situation, we force RHOBEG and RHOEND to equal when the difference is tiny.
          rhoend = rhobeg
      end if

      if (rhobeg <= ZERO .or. is_nan(rhobeg) .or. is_inf(rhobeg))then
          rhobeg = max(TEN*RHOEND, RHOBEG_DFT)
          print '(/1X, 1A, 1PD15.6, 1A)', solver // ': invalid RHOBEG; i&
     &t should be a positive number; it is set to ', rhobeg, '.'
      end if

      if (rhoend < ZERO .or. rhobeg < rhoend .or. is_nan(rhoend) .or. is&
     &_inf(rhoend)) then
          rhoend = max(EPS, min(TENTH*rhobeg, RHOEND_DFT))
          print '(/1X, 1A, 1PD15.6, 1A)', solver // ': invalid RHOEND; i&
     &t should be a positive number and RHOEND <= RHOBEG; ' // 'it is se&
     &t to ', rhoend, '.'
      end if

      end subroutine


      end module preproc_mod