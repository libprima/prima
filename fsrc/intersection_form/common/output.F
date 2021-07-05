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
! on 04-Jul-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! output_mod is a module providing some subroutines concerning output to
! terminla/files. Note that these output operations are sequential in
! nature. In case parallelisum is desirable (especially during
! initializaton), the subroutines may have to be modified or disabled.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and papers.
!
! Last Modified: Sunday, May 23, 2021 AM11:08:31


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
      character(len=*), intent(in) :: solver

! Local variables
      integer :: ios
      ! Should be an integer of default kind
      character(len=100) :: fout
      character(len=3) :: fstat
      character(len=MSSGLEN) :: mssg
      logical :: fexist


      if (info == FTARGET_ACHIEVED) then
          mssg = 'the target function value is achieved.'
      else if (info == MAXFUN_REACHED) then
          mssg = 'the objective function has been evaluated MAXFUN times&
     &.'
      else if (info == SMALL_TR_RADIUS) then
          mssg = 'the trust region radius reaches its lower bound.'
      else if (info == TRSUBP_FAILED) then
          mssg = 'a trust region step has failed to reduce the quadratic&
     & model.'
      else if (info == NAN_X) then
          mssg = 'NaN occurs in x.'
      else if (info == NAN_INF_F) then
          mssg = 'the objective function returns NaN or +INFINITY.'
      else if (info == NAN_MODEL) then
          mssg = 'NaN occurs in the models.'
      end if

      if (iprint >= 1) then
          if (iprint >= 3) then
              print '(1X)'
          end if
          print '(/4A)', 'Return from ', solver, ' because ', trim(mssg)
          print '(1A, 6X, 1A, I7)', 'At the return from '//solver, 'Numb&
     &er of function evaluations = ', nf
          print '(1A, 1PD23.15, 6X, 1A, /(1P, 5D15.6))', 'Least function&
     & value = ', f, 'The corresponding X is:', x
          print '(1X)'
      end if

      if (iprint <= -1) then
          fout = solver//'_output.txt'
          inquire (file=trim(fout), exist=fexist)
          if (fexist) then
              fstat = 'old'
          else
              fstat = 'new'
          end if
          open (unit=OUTUNIT, file=trim(fout), status=fstat, position='a&
     &ppend', iostat=ios, action='write')
          if (ios /= 0) then
              print '(1A)', 'Fail to open file '//trim(fout)//'!'
          else
              if (iprint <= -3) then
                  write (OUTUNIT, '(1X)')
              end if
              write (OUTUNIT, '(/4A)') 'Return from ', solver, ' because&
     & ', trim(mssg)
              write (OUTUNIT, '(1A, 6X, 1A, I7)') 'At the return from '/&
     &/solver, 'Number of function evaluations = ', nf
              write (OUTUNIT, '(1A, 1PD23.15, 6X, 1A, /(1P, 5D15.6))') '&
     &Least function value = ', f, 'The corresponding X is:', x
              close (OUTUNIT)
          end if
!print '(/1A /)', 'The output is printed to ' // trim(fout) // '.'
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
      character(len=*), intent(in) :: solver

! Local variables
      integer :: ios
      ! Should be an integer of default kind
      character(len=100) :: fout
      character(len=3) :: fstat
      logical :: fexist


      if (iprint >= 2) then
          if (iprint >= 3) then
              print '(1X)'
          end if
          print '(/1A, 1PD11.4, 6X, A, I7)', 'New RHO = ', rho, 'Number &
     &of function evaluations = ', nf
          print '(1A, 1PD23.15, 6X, 1A, /(1P, 5D15.6))', 'Least function&
     & value = ', f, 'The corresponding X is:', x
      end if

      if (iprint <= -2) then
          fout = solver//'_output.txt'
          inquire (file=trim(fout), exist=fexist)
          if (fexist) then
              fstat = 'old'
          else
              fstat = 'new'
          end if
          open (unit=OUTUNIT, file=trim(fout), status=fstat, position='a&
     &ppend', iostat=ios, action='write')
          if (ios /= 0) then
              print '(1A)', 'Fail to open file '//trim(fout)//'!'
          else
              if (iprint <= -3) then
                  write (OUTUNIT, '(1X)')
              end if
              write (OUTUNIT, '(/1A, 1PD11.4, 6X, A, I7)') 'New RHO = ',&
     & rho, 'Number of function evaluations = ', nf
              write (OUTUNIT, '(1A, 1PD23.15, 6X, 1A, /(1P, 5D15.6))') '&
     &Least function value = ', f, 'The corresponding X is:', x
              close (OUTUNIT)
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
      character(len=*), intent(in) :: solver

! Local variables
      integer :: ios
      ! Should be an integer of default kind
      character(len=100) :: fout
      character(len=3) :: fstat
      logical :: fexist

!Function number     4    F =  1.1853827160D-01    The corresponding X is:

      if (iprint >= 3) then
          print '(/1A, I7, 4X, 1A, 1PD18.10, 4X, 1A, /(1P, 5D15.6))', 'F&
     &unction number', nf, 'F = ', f, 'The corresponding X is:', x
      end if

      if (iprint <= -3) then
          fout = solver//'_output.txt'
          inquire (file=trim(fout), exist=fexist)
          if (fexist) then
              fstat = 'old'
          else
              fstat = 'new'
          end if
          open (unit=OUTUNIT, file=trim(fout), status=fstat, position='a&
     &ppend', iostat=ios, action='write')
          if (ios /= 0) then
              print '(1A)', 'Fail to open file '//trim(fout)//'!'
          else
              write (OUTUNIT, '(/1A, I7, 4X, 1A, 1PD18.10, 4X, 1A, /(1P,&
     & 5D15.6))') 'Function number', nf, 'F = ', f, 'The corresponding X&
     & is:', x
              close (OUTUNIT)
          end if
      end if

      end subroutine fmssg


      end module output_mod