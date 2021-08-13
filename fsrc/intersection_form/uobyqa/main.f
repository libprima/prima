!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of main.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 14-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==aa0001.f90  processed by SPAG 7.50RE at 00:28 on 26 May 2021
!
!     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6,8.
!
            implicit none
!*--AA00018
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            integer :: i, iprint, maxfun, n, info
            real*8(R8KIND) :: rhobeg, rhoend, f, ftarget
            real*8(R8KIND), dimension(10000) :: w
            real*8(R8KIND), dimension(10) :: x
!*++
!*++ End of declarations rewritten by SPAG
!*++
            iprint = 2
            maxfun = 5000
            rhoend = 1.0D-8
            do n = 2, 8, 2
                do i = 1, n
                    x(i) = DFLOAT(i) / DFLOAT(n + 1)
                end do
                rhobeg = 0.2D0 * x(1)
                print 99001, n
99001 format(//5X, '******************'/5X, 'Results with N =', I2, /5X,&
     & '******************')
                call UOBYQA(n, x, rhobeg, rhoend, iprint, maxfun, w, f, &
     &info, ftarget)
            end do
            end program AA0001