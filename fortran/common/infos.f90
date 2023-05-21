module infos_mod
!--------------------------------------------------------------------------------------------------!
! This is a module defining exit flags.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Sunday, May 21, 2023 PM03:00:41
!--------------------------------------------------------------------------------------------------!

use, non_intrinsic :: consts_mod, only : IK

implicit none
private
public :: INFO_DFT
public :: SMALL_TR_RADIUS
public :: FTARGET_ACHIEVED
public :: TRSUBP_FAILED
public :: MAXFUN_REACHED
public :: MAXTR_REACHED
public :: NAN_INF_X
public :: NAN_INF_F
public :: NAN_INF_MODEL
public :: DAMAGING_ROUNDING
public :: NO_SPACE_BETWEEN_BOUNDS
public :: ZERO_LINEAR_CONSTRAINT
public :: INVALID_INPUT
public :: ASSERTION_FAILS
public :: VALIDATION_FAILS
public :: MEMORY_ALLOCATION_FAILS


integer(IK), parameter :: INFO_DFT = 0
integer(IK), parameter :: SMALL_TR_RADIUS = 0
integer(IK), parameter :: FTARGET_ACHIEVED = 1
integer(IK), parameter :: TRSUBP_FAILED = 2
integer(IK), parameter :: MAXFUN_REACHED = 3
integer(IK), parameter :: MAXTR_REACHED = 20
integer(IK), parameter :: NAN_INF_X = -1
integer(IK), parameter :: NAN_INF_F = -2
integer(IK), parameter :: NAN_INF_MODEL = -3
integer(IK), parameter :: NO_SPACE_BETWEEN_BOUNDS = 6
integer(IK), parameter :: DAMAGING_ROUNDING = 7
integer(IK), parameter :: ZERO_LINEAR_CONSTRAINT = 8

! Stop-codes.
! The following codes are used by ERROR STOP as stop-codes, which should be default integers.
integer, parameter :: INVALID_INPUT = 100
integer, parameter :: ASSERTION_FAILS = 101
integer, parameter :: VALIDATION_FAILS = 102
integer, parameter :: MEMORY_ALLOCATION_FAILS = 103


end module infos_mod
