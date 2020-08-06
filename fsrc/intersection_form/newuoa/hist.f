!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of hist.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 07-Aug-2020.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      module hist

      updatehist(ihist, )

! Inputs
      integer(IK) :: ihist
      integer(IK) :: nf
      real(RP) :: f
      real(RP) :: x(:)


! In-outputs
      real(RP) :: fhist(:)
      real(RP) :: xhist(:, :)

      k = mod(nf - 1, size(fhist)) + 1
      fhist(k) = f
      xhist(:, k) = x

      end if


      end module hist