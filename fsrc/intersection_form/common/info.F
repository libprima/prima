!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of info.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 23-May-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! INFO_MOD is a module defining exit flags.
!
! Coded by Zaikun Zhang in July 2020.
!
! Last Modified: Sunday, May 23, 2021 AM11:08:02


      module info_mod


      implicit none
      private
      public :: INVALID_INPUT
      public :: SMALL_TR_RADIUS
      public :: FTARGET_ACHIEVED
      public :: TRSUBP_FAILED
      public :: MAXFUN_REACHED
      public :: NAN_X
      public :: NAN_INF_F
      public :: NAN_MODEL

      integer, parameter :: INVALID_INPUT = -100
      integer, parameter :: SMALL_TR_RADIUS = 0
      integer, parameter :: FTARGET_ACHIEVED = 1
      integer, parameter :: TRSUBP_FAILED = 2
      integer, parameter :: MAXFUN_REACHED = 3
      integer, parameter :: NAN_X = -1
      integer, parameter :: NAN_INF_F = -2
      integer, parameter :: NAN_MODEL = -3


      end module info_mod