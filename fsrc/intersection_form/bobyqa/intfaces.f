!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of intfaces.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 25-May-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==intfaces.f90  created by SPAG 7.50RE at 17:55 on 25 May 2021
      module S_ALTMOV
      interface
          subroutine ALTMOV(N, Npt, Xpt, Xopt, Bmat, Zmat, Ndim, Sl, Su,&
     & Kopt, Knew, Adelt, Xnew, Xalt, Alpha, Cauchy, Glag, Hcol, W)
          use F77KINDS
          implicit none
          integer, intent(IN) :: N
          integer, intent(IN) :: Npt
          real*8, intent(IN), dimension(Npt, *) :: Xpt
          real*8, intent(IN), dimension(*) :: Xopt
          real*8, intent(IN), dimension(Ndim, *) :: Bmat
          real*8, intent(IN), dimension(Npt, *) :: Zmat
          integer, intent(IN) :: Ndim
          real*8, intent(IN), dimension(*) :: Sl
          real*8, intent(IN), dimension(*) :: Su
          integer, intent(IN) :: Kopt
          integer, intent(IN) :: Knew
          real*8, intent(IN) :: Adelt
          real*8, intent(OUT), dimension(*) :: Xnew
          real*8, intent(INOUT), dimension(*) :: Xalt
          real*8, intent(INOUT) :: Alpha
          real*8, intent(INOUT) :: Cauchy
          real*8, intent(INOUT), dimension(*) :: Glag
          real*8, intent(INOUT), dimension(*) :: Hcol
          real*8, intent(INOUT), dimension(*) :: W
          end subroutine ALTMOV
      end interface
      end module S_ALTMOV
!*==intfaces.f90  created by SPAG 7.50RE at 17:55 on 25 May 2021
      module S_BOBYQA
      interface
          subroutine BOBYQA(N, Npt, X, Xl, Xu, Rhobeg, Rhoend, Iprint, M&
     &axfun, W, F, Info, Ftarget)
          use F77KINDS
          implicit none
          integer :: N
          integer :: Npt
          real*8, intent(INOUT), dimension(*) :: X
          real*8, dimension(*) :: Xl
          real*8, dimension(*) :: Xu
          real*8, intent(INOUT) :: Rhobeg
          real*8, intent(INOUT) :: Rhoend
          integer :: Iprint
          integer :: Maxfun
          real*8, intent(INOUT), dimension(*) :: W
          real*8 :: F
          integer :: Info
          real*8 :: Ftarget
          end subroutine BOBYQA
      end interface
      end module S_BOBYQA
!*==intfaces.f90  created by SPAG 7.50RE at 17:55 on 25 May 2021
      module S_BOBYQB
      interface
          subroutine BOBYQB(N, Npt, X, Xl, Xu, Rhobeg, Rhoend, Iprint, M&
     &axfun, Xbase, Xpt, Fval, Xopt, Gopt, Hq, Pq, Bmat, Zmat, Ndim, Sl,&
     & Su, Xnew, Xalt, D, Vlag, W, F, Info, Ftarget)
          use F77KINDS
          implicit none
          integer :: N
          integer :: Npt
          real*8, intent(INOUT), dimension(*) :: X
          real*8, dimension(*) :: Xl
          real*8, dimension(*) :: Xu
          real*8 :: Rhobeg
          real*8, intent(IN) :: Rhoend
          integer :: Iprint
          integer :: Maxfun
          real*8, intent(INOUT), dimension(*) :: Xbase
          real*8, intent(INOUT), dimension(Npt, *) :: Xpt
          real*8, intent(INOUT), dimension(*) :: Fval
          real*8, intent(INOUT), dimension(*) :: Xopt
          real*8, intent(INOUT), dimension(*) :: Gopt
          real*8, intent(INOUT), dimension(*) :: Hq
          real*8, intent(INOUT), dimension(*) :: Pq
          real*8, intent(INOUT), dimension(Ndim, *) :: Bmat
          real*8, dimension(Npt, *) :: Zmat
          integer :: Ndim
          real*8, intent(INOUT), dimension(*) :: Sl
          real*8, intent(INOUT), dimension(*) :: Su
          real*8, intent(INOUT), dimension(*) :: Xnew
          real*8, dimension(*) :: Xalt
          real*8, intent(INOUT), dimension(*) :: D
          real*8, intent(INOUT), dimension(*) :: Vlag
          real*8, intent(INOUT), dimension(*) :: W
          real*8, intent(INOUT) :: F
          integer, intent(OUT) :: Info
          real*8 :: Ftarget
          end subroutine BOBYQB
      end interface
      end module S_BOBYQB
!*==intfaces.f90  created by SPAG 7.50RE at 17:55 on 25 May 2021
      module S_CALFUN
      interface
          subroutine CALFUN(N, X, F)
          use F77KINDS
          implicit none
          integer, intent(IN) :: N
          real*8, intent(IN), dimension(*) :: X
          real*8, intent(INOUT) :: F
          end subroutine CALFUN
      end interface
      end module S_CALFUN
!*==intfaces.f90  created by SPAG 7.50RE at 17:55 on 25 May 2021
      module S_PRELIM
      interface
          subroutine PRELIM(N, Npt, X, Xl, Xu, Rhobeg, Iprint, Maxfun, X&
     &base, Xpt, Fval, Gopt, Hq, Pq, Bmat, Zmat, Ndim, Sl, Su, Nf, Kopt,&
     & F, Ftarget)
          use F77KINDS
          implicit none
          integer :: N
          integer, intent(IN) :: Npt
          real*8, intent(INOUT), dimension(*) :: X
          real*8, intent(IN), dimension(*) :: Xl
          real*8, intent(IN), dimension(*) :: Xu
          real*8, intent(IN) :: Rhobeg
          integer, intent(IN) :: Iprint
          integer, intent(IN) :: Maxfun
          real*8, intent(INOUT), dimension(*) :: Xbase
          real*8, intent(INOUT), dimension(Npt, *) :: Xpt
          real*8, intent(INOUT), dimension(*) :: Fval
          real*8, intent(INOUT), dimension(*) :: Gopt
          real*8, intent(OUT), dimension(*) :: Hq
          real*8, intent(OUT), dimension(*) :: Pq
          real*8, intent(INOUT), dimension(Ndim, *) :: Bmat
          real*8, intent(INOUT), dimension(Npt, *) :: Zmat
          integer, intent(IN) :: Ndim
          real*8, intent(IN), dimension(*) :: Sl
          real*8, intent(IN), dimension(*) :: Su
          integer, intent(INOUT) :: Nf
          integer, intent(INOUT) :: Kopt
          real*8 :: F
          real*8, intent(IN) :: Ftarget
          end subroutine PRELIM
      end interface
      end module S_PRELIM
!*==intfaces.f90  created by SPAG 7.50RE at 17:55 on 25 May 2021
      module S_RESCUE
      interface
          subroutine RESCUE(N, Npt, Xl, Xu, Iprint, Maxfun, Xbase, Xpt, &
     &Fval, Xopt, Gopt, Hq, Pq, Bmat, Zmat, Ndim, Sl, Su, Nf, Delta, Kop&
     &t, Vlag, Ptsaux, Ptsid, W, F, Ftarget)
          use F77KINDS
          implicit none
          integer :: N
          integer :: Npt
          real*8, intent(IN), dimension(*) :: Xl
          real*8, intent(IN), dimension(*) :: Xu
          integer, intent(IN) :: Iprint
          integer, intent(IN) :: Maxfun
          real*8, intent(INOUT), dimension(*) :: Xbase
          real*8, intent(INOUT), dimension(Npt, *) :: Xpt
          real*8, intent(INOUT), dimension(*) :: Fval
          real*8, intent(INOUT), dimension(*) :: Xopt
          real*8, intent(INOUT), dimension(*) :: Gopt
          real*8, intent(INOUT), dimension(*) :: Hq
          real*8, intent(INOUT), dimension(*) :: Pq
          real*8, intent(INOUT), dimension(Ndim, *) :: Bmat
          real*8, intent(INOUT), dimension(Npt, *) :: Zmat
          integer :: Ndim
          real*8, intent(INOUT), dimension(*) :: Sl
          real*8, intent(INOUT), dimension(*) :: Su
          integer, intent(INOUT) :: Nf
          real*8, intent(IN) :: Delta
          integer, intent(INOUT) :: Kopt
          real*8, intent(INOUT), dimension(*) :: Vlag
          real*8, intent(INOUT), dimension(2, *) :: Ptsaux
          real*8, intent(INOUT), dimension(*) :: Ptsid
          real*8, intent(INOUT), dimension(*) :: W
          real*8 :: F
          real*8, intent(IN) :: Ftarget
          end subroutine RESCUE
      end interface
      end module S_RESCUE
!*==intfaces.f90  created by SPAG 7.50RE at 17:55 on 25 May 2021
      module S_TRSBOX
      interface
          subroutine TRSBOX(N, Npt, Xpt, Xopt, Gopt, Hq, Pq, Sl, Su, Del&
     &ta, Xnew, D, Gnew, Xbdi, S, Hs, Hred, Dsq, Crvmin)
          use F77KINDS
          implicit none
          integer, intent(IN) :: N
          integer, intent(IN) :: Npt
          real*8, intent(IN), dimension(Npt, *) :: Xpt
          real*8, intent(IN), dimension(*) :: Xopt
          real*8, intent(IN), dimension(*) :: Gopt
          real*8, intent(IN), dimension(*) :: Hq
          real*8, intent(IN), dimension(*) :: Pq
          real*8, intent(IN), dimension(*) :: Sl
          real*8, intent(IN), dimension(*) :: Su
          real*8, intent(IN) :: Delta
          real*8, intent(INOUT), dimension(*) :: Xnew
          real*8, intent(INOUT), dimension(*) :: D
          real*8, intent(INOUT), dimension(*) :: Gnew
          real*8, intent(INOUT), dimension(*) :: Xbdi
          real*8, intent(INOUT), dimension(*) :: S
          real*8, intent(INOUT), dimension(*) :: Hs
          real*8, intent(INOUT), dimension(*) :: Hred
          real*8, intent(INOUT) :: Dsq
          real*8, intent(INOUT) :: Crvmin
          end subroutine TRSBOX
      end interface
      end module S_TRSBOX
!*==intfaces.f90  created by SPAG 7.50RE at 17:55 on 25 May 2021
      module S_UPDATE
      interface
          subroutine UPDATE(N, Npt, Bmat, Zmat, Ndim, Vlag, Beta, Denom,&
     & Knew, W)
          use F77KINDS
          implicit none
          integer, intent(IN) :: N
          integer, intent(IN) :: Npt
          real*8, intent(INOUT), dimension(Ndim, *) :: Bmat
          real*8, intent(INOUT), dimension(Npt, *) :: Zmat
          integer, intent(IN) :: Ndim
          real*8, intent(INOUT), dimension(*) :: Vlag
          real*8, intent(IN) :: Beta
          real*8, intent(IN) :: Denom
          integer, intent(IN) :: Knew
          real*8, intent(INOUT), dimension(*) :: W
          end subroutine UPDATE
      end interface
      end module S_UPDATE