!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of dirty_temporary_mod4powell.f90.
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


      module dirty_temporary_mod4powell_mod
!--------------------------------------------------------------------------------------------------!
! This is a DIRTY and TEMPORARY module to make some constants and subroutines available in Powell's
! original code. It is NEVER used in the modernized code.
!--------------------------------------------------------------------------------------------------!
      use consts_mod, only : ZERO, ONE, TWO, HALF, TEN, TENTH, QUART, PI
      use consts_mod, only : HUGENUM, HUGEFUN, HUGECON
      use consts_mod, only : IK, RP
      use linalg_mod, only : matprod, inprod, norm, calquad, inprod, ism&
     &inor, planerot, eye, hypotenuse, project, inv
      use linalg_mod, only : matmul => matprod, dot_product => inprod
      use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED, TRSUBP_FAIL&
     &ED, SMALL_TR_RADIUS, NAN_INF_X, NAN_INF_F
      use infnan_mod, only : is_nan, is_posinf, is_inf, is_posinf, is_ne&
     &ginf
      use ratio_mod, only : redrat
      use redrho_mod, only : redrho
!use debug_mod, only : errstop, assert
!use output_mod, only : retmsg, rhomsg, fmsg

      end module dirty_temporary_mod4powell_mod