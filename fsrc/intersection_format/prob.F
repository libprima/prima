! This is the intersection-format version of /home/zaikunzhang/Work/neupdfo/fsrc/intersection_format/prob.F.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection format, each continued line has an ampersand at
! column 73, and each continuation line has an ampersand at column 6.
! A Fortran file in such a format can be compiled both as a fixed-format
! file and as a free-format file.
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 03-Aug-2020.


! PROB_MOD is a module specifying the interfaces for CALFUN and CALCFC.
! CALFUN evaluates the objective function for unconstrained, bound
! constrained, and linearly constrained problems; CALCFC evaluates the
! objective function and constraint for nonlinearly constrained prolems.
!
! Users must provide the implementation of CALFUN or CALCFC.

            module prob_mod

            implicit none
            private

            public :: calfun, calcfc

            interface
                subroutine calfun(x, f)
                use consts_mod, only : RP
                implicit none
                real(RP), intent(in) :: x(:)
                real(RP), intent(out) :: f
                end subroutine calfun

                subroutine calcfc(x, f, con)
                use consts_mod, only : RP
                implicit none
                real(RP), intent(in) :: x(:)
                real(RP), intent(out) :: f
                real(RP), intent(out) :: con(:)
                end subroutine calcfc
            end interface

            end module prob_mod