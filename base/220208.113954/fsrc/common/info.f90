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
