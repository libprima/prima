module dirty_temporary_mod4powell_mod
!--------------------------------------------------------------------------------------------------!
! This is a DIRTY and TEMPORARY module to make some constants and subroutines available in Powell's
! original code. It is NEVER used in the modernized code.
!--------------------------------------------------------------------------------------------------!
use consts_mod, only : ZERO, ONE, TWO, HALF, TEN, TENTH, QUART, PI
use consts_mod, only : HUGENUM, HUGEFUN, HUGECON
use consts_mod, only : IK, RP
use linalg_mod, only : matprod, inprod, norm, calquad, inprod, isminor, planerot, eye, hypotenuse, project, inv
use linalg_mod, only : matmul => matprod, dot_product => inprod
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAILED, SMALL_TR_RADIUS, NAN_INF_X, NAN_INF_F
use infnan_mod, only : is_nan, is_posinf, is_inf, is_posinf, is_neginf
!use debug_mod, only : errstop, assert
!use output_mod, only : retmsg, rhomsg, fmsg

end module dirty_temporary_mod4powell_mod
