!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of prob.F90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 08-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! PROB_MOD is a module specifying the interfaces for CALFUN and CALCFC.
! CALFUN evaluates the objective function for unconstrained, bound
! constrained, and linearly constrained problems; CALCFC evaluates the
! objective function and constraint for nonlinearly constrained prolems.
!
! Coded by Zaikun Zhang in July 2020.
!
! Users must provide the implementation of CALFUN or CALCFC.

      module prob_mod

      implicit none
      private

      public :: funeval, fceval

      abstract interface
          subroutine funeval(x, f)
          use consts_mod, only : RP
          implicit none
          real(RP), intent(in) :: x(:)
          real(RP), intent(out) :: f
          end subroutine funeval
      end interface


      abstract interface
          subroutine fceval(x, f, con)
          use consts_mod, only : RP
          implicit none
          real(RP), intent(in) :: x(:)
          real(RP), intent(out) :: f
          real(RP), intent(out) :: con(:)
          end subroutine fceval
      end interface

      end module prob_mod