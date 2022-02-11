!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of pintrf.f90.
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


      module pintrf_mod
!--------------------------------------------------------------------------------------------------!
! This is a module specifying the abstract interfaces OBJ and OBJCON. OBJ evaluates the objective
! function for unconstrained, bound constrained, and linearly constrained problems; OBJCON evaluates
! the objective and constraint functions for nonlinearly constrained problems.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: July 2020.
!
! Last Modified: Monday, February 07, 2022 AM12:28:54
!--------------------------------------------------------------------------------------------------!

!!!!!! Users must provide the implementation of OBJ or OBJCON. !!!!!!

      implicit none
      private
      public :: OBJ, OBJCON


      abstract interface

          subroutine OBJ(x, f)
          use consts_mod, only : RP
          implicit none
          real(RP), intent(in) :: x(:)
          real(RP), intent(out) :: f
          end subroutine OBJ


          subroutine OBJCON(x, f, constr)
          use consts_mod, only : RP
          implicit none
          real(RP), intent(in) :: x(:)
          real(RP), intent(out) :: f
          real(RP), intent(out) :: constr(:)
          end subroutine OBJCON

      end interface


      end module pintrf_mod