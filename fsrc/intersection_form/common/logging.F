!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of logging.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 06-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! logging_mod is a module providing some subroutines for logging.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and papers.
!
! Last Modified: Friday, July 02, 2021 AM11:39:27


      module logging_mod

      implicit none
      private
      public :: logging


      contains


      subroutine logging(logfile, srname, lnnum, nf, f, x, con, conv, ms&
     &sg)
      use consts_mod, only : RP, IK
      implicit none

! Inputs
      integer(IK), intent(in) :: lnnum
      integer(IK), intent(in) :: nf
      real(RP), intent(in) :: f
      real(RP), optional, intent(in) :: x(:)
      real(RP), optional, intent(in) :: con(:)
      real(RP), optional, intent(in) :: conv
      character(len=*), intent(in) :: logfile
      character(len=*), intent(in) :: srname
      character(len=*), optional, intent(in) :: mssg

! Local variables
      integer(IK), parameter :: LOGUNIT = 11
      integer :: ios
      ! Should be an integer of default kind
      character(len=3) :: fstat
      logical :: fexist


      inquire (file=trim(logfile), exist=fexist)
      if (fexist) then
          fstat = 'old'
      else
          fstat = 'new'
      end if
      open (unit=LOGUNIT, file=trim(logfile), status=fstat, position='ap&
     &pend', iostat=ios, action='write')
      if (ios /= 0) then
          print '(1A)', 'Fail to open file '//trim(logfile)//'!'
      else
          write (LOGUNIT, '(/1A, I5, 4X, 1A, I7, 4X, 1A, 1PD23.15)') 'In&
     & '//trim(srname)//', Line number', lnnum, 'Function number', nf, '&
     &F = ', f
          if (present(con)) then
              write (LOGUNIT, '(1A, 1PD23.15)') 'The constraint violatio&
     &n is:', conv
          end if
          if (present(mssg)) then
              write (LOGUNIT, '(1A)') 'Message: '//mssg
          end if
          if (present(x)) then
              write (LOGUNIT, '(1A, /(1P, 5D23.15))') 'The corresponding&
     & X is:', x
          end if
          if (present(con)) then
              write (LOGUNIT, '(1A, /(1P, 5D23.15))') 'The constraint is&
     &:', con
          end if
          close (LOGUNIT)
      end if

      end subroutine logging


      end module logging_mod