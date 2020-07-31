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
