!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of info.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun ZHANG (www.zhangzk.net)
! on 11-Feb-2022.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! INFO_MOD is a module defining exit flags.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Sunday, September 26, 2021 AM08:29:06


      module info_mod

      use, non_intrinsic :: consts_mod, only : IK

      implicit none
      private
      public :: INFO_DFT
      public :: INVALID_INPUT
      public :: SMALL_TR_RADIUS
      public :: FTARGET_ACHIEVED
      public :: TRSUBP_FAILED
      public :: MAXFUN_REACHED
      public :: MAXTR_REACHED
      public :: NAN_INF_X
      public :: NAN_INF_F
      public :: NAN_MODEL
      public :: DAMAGING_ROUNDING


      integer(IK), parameter :: INFO_DFT = 0_IK
      integer(IK), parameter :: INVALID_INPUT = -100_IK
      integer(IK), parameter :: SMALL_TR_RADIUS = 0_IK
      integer(IK), parameter :: FTARGET_ACHIEVED = 1_IK
      integer(IK), parameter :: TRSUBP_FAILED = 2_IK
      integer(IK), parameter :: MAXFUN_REACHED = 3_IK
      integer(IK), parameter :: MAXTR_REACHED = 20_IK
      integer(IK), parameter :: NAN_INF_X = -1_IK
      integer(IK), parameter :: NAN_INF_F = -2_IK
      integer(IK), parameter :: NAN_MODEL = -3_IK
      integer(IK), parameter :: DAMAGING_ROUNDING = 7_IK


      end module info_mod