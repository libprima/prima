!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of output.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 06-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! output_mod is a module providing some subroutines concerning output to
! terminla/files. Note that these output operations are sequential in
! nature. In case parallelisum is desirable (especially during
! initializaton), the subroutines may have to be modified or disabled.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code
! and the NEWUOA paper.

      module output_mod

      implicit none
      private
      public :: retmssg, rhomssg, fmssg


      contains


      subroutine retmssg(info, iprint, nf, f, x, solver)
      use consts_mod, only : RP, IK, MSSGLEN, OUTUNIT
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED
      use info_mod, only : SMALL_TR_RADIUS, TRSUBP_FAILED
      use info_mod, only : NAN_X, NAN_INF_F, NAN_MODEL
      implicit none

! Inputs
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: info
      integer(IK), intent(in) :: nf
      real(RP), intent(in) :: f
      real(RP), intent(in) :: x(:)
      character(len = *), intent(in) :: solver

! Intermediate variables
      integer :: ios ! Should be an integer of default kind
      character(len = 100) :: fout
      character(len = 3) :: fstat
      character(len = MSSGLEN) :: mssg
      logical :: fexist


      if (info == FTARGET_ACHIEVED) then
          mssg = 'the target function value is achieved.'
      else if (info == MAXFUN_REACHED) then
          mssg = 'the objective function has been evaluated MAXFUN times&
     &.'
      else if (info == SMALL_TR_RADIUS) then
          mssg = 'the lower bound for the trust region radius is reached&
     &.'
      else if (info == TRSUBP_FAILED) then
          mssg = 'a trust region step has failed to reduce the quadratic&
     &model.'
      else if (info == NAN_X) then
          mssg = 'NaN occurs in x.'
      else if (info == NAN_INF_F) then
          mssg = 'the objective function returns NaN or +INFINITY.'
      else if (info == NAN_MODEL) then
          mssg = 'NaN occurs in the models.'
      end if

      if (iprint >= 1) then
          print '(/1X, 4A)', 'Return from ', solver, ' because ', trim(m&
     &ssg)
          print '(1X, 1A, I7)', 'Number of function evaluations = ', nf
          print '(1X, 1A, 1PD23.15)', 'Least function value = ', f
          print '(1X, 1A, /(1X, 1P, 5D15.6))', 'The corresponding X is',&
     &x
      end if

      if (iprint >= 4) then
          fout = solver // '.output'
          inquire(file = trim(fout), exist = fexist)
          if (fexist) then
              fstat = 'old'
          else
              fstat = 'new'
          end if
          open(unit = OUTUNIT, file = trim(fout), status = fstat, positi&
     &on = 'append', iostat = ios, action = 'write')
          if (ios /= 0) then
              print '(1X, 1A)', 'Fail to open file ' // fout // '!'
          else
              write(OUTUNIT, '(/1X, 4A)') 'Return from ', solver, ' beca&
     &use ', trim(mssg)
              write(OUTUNIT, '(1X, 1A, I7)') 'Number of function evaluat&
     &ions = ', nf
              write(OUTUNIT, '(1X, 1A, 1PD23.15)') 'Least function value&
     &= ', f
              write(OUTUNIT, '(1X, 1A, /(1X, 1P, 5D15.6))') 'The corresp&
     &onding X is', x
              close(OUTUNIT)
          end if
      end if

      end subroutine retmssg


      subroutine rhomssg(iprint, nf, f, rho, x, solver)
      use consts_mod, only : RP, IK, OUTUNIT
      implicit none

! Inputs
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: nf
      real(RP), intent(in) :: f
      real(RP), intent(in) :: rho
      real(RP), intent(in) :: x(:)
      character(len = *), intent(in) :: solver

! Intermediate variables
      integer :: ios ! Should be an integer of default kind
      character(len = 100) :: fout
      character(len = 3) :: fstat
      logical :: fexist


      if (iprint >= 2) then
          print '(/1X, 1A, 1PD11.4)', 'New RHO = ', rho
          print '(1X, 1A, I7)', 'Number of function evaluations = ', nf
          print '(1X, 1A, 1PD23.15)', 'Least function value = ', f
          print '(1X, 1A, /(1X, 1P, 5D15.6))', 'The corresponding X is',&
     &x
      end if

      if (iprint >= 4) then
          fout = solver // '.output'
          inquire(file = trim(fout), exist = fexist)
          if (fexist) then
              fstat = 'old'
          else
              fstat = 'new'
          end if
          open(unit = OUTUNIT, file = trim(fout), status = fstat, positi&
     &on = 'append', iostat = ios, action = 'write')
          if (ios /= 0) then
              print '(1X, 1A)', 'Fail to open file ' // fout // '!'
          else
              write(OUTUNIT, '(/1X, 1A, 1PD11.4)') 'New RHO = ', rho
              write(OUTUNIT, '(1X, 1A, I7)') 'Number of function evaluat&
     &ions = ', nf
              write(OUTUNIT, '(1X, 1A, 1PD23.15)') 'Least function value&
     &= ', f
              write(OUTUNIT, '(1X, 1A, /(1X, 1P, 5D15.6))') 'The corresp&
     &onding X is', x
              close(OUTUNIT)
          end if
      end if

      end subroutine rhomssg


      subroutine fmssg(iprint, nf, f, x, solver)
      use consts_mod, only : RP, IK, OUTUNIT
      implicit none

! Inputs
      integer(IK), intent(in) :: iprint
      integer(IK), intent(in) :: nf
      real(RP), intent(in) :: f
      real(RP), intent(in) :: x(:)
      character(len = *), intent(in) :: solver

! Intermediate variables
      integer :: ios ! Should be an integer of default kind
      character(len = 100) :: fout
      character(len = 3) :: fstat
      logical :: fexist


      if (iprint >= 3) then
         print '(/1X, 1A, I7)', 'Function number', nf
         print '(1X, 1A, 1PD23.15)', 'Function value = ', f
         print '(1X, 1A, /(1X, 1P, 5D15.6))', 'The corresponding X is:',&
     &x
      end if

      if (iprint >= 4) then
          fout = solver // '.output'
          inquire(file = trim(fout), exist = fexist)
          if (fexist) then
              fstat = 'old'
          else
              fstat = 'new'
          end if
          open(unit = OUTUNIT, file = trim(fout), status = fstat, positi&
     &on = 'append', iostat = ios, action = 'write')
          if (ios /= 0) then
              print '(1X, 1A)', 'Fail to open file ' // fout // '!'
          else
              write(OUTUNIT, '(/1X, 1A, I7)') 'Function number', nf
              write(OUTUNIT, '(1X, 1A, 1PD23.15)') 'Function value = ', &
     &f
              write(OUTUNIT, '(1X, 1A, /(1X, 1P, 5D15.6))') 'The corresp&
     &onding X is:', x
              close(OUTUNIT)
          end if
      end if

      end subroutine fmssg


      end module output_mod