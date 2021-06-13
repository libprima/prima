!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of prelim.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 14-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==prelim.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
            subroutine PRELIM(N, Npt, X, Xl, Xu, Rhobeg, Iprint, Maxfun,&
     & Xbase, Xpt, Fval, Gopt, Hq, Pq, Bmat, Zmat, Ndim, Sl, Su, Nf, Kop&
     &t, F, Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--PRELIM13
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
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
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            real*8 :: almost_infinity, diff, fbeg, half, one, recip, rho&
     &sq, stepa, stepb, temp, two, zero
            integer :: i, ih, ipt, itemp, j, jpt, k, nfm, nfx, np
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonzero, BOBYQB will change it to its usual value later.
!     NF is maintaned as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!
!     Set some constants.
!
            half = 0.5D0
            one = 1.0D0
            two = 2.0D0
            zero = 0.0D0
            rhosq = Rhobeg * Rhobeg
            recip = one / rhosq
            np = N + 1
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            almost_infinity = huge(0.0D0) / 2.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
            do j = 1, N
                Xbase(j) = X(j)
                do k = 1, Npt
                    Xpt(k, j) = zero
                end do
                do i = 1, Ndim
                    Bmat(i, j) = zero
                end do
            end do
            do ih = 1, (N * np) / 2
                Hq(ih) = zero
            end do
            do k = 1, Npt
                Pq(k) = zero
                do j = 1, Npt - np
                    Zmat(k, j) = zero
                end do
            end do
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
            Nf = 0
            do
                nfm = Nf
                nfx = Nf - N
                Nf = Nf + 1
                if (nfm > 2 * N) then
                    itemp = (nfm - np) / N
                    jpt = nfm - itemp * N - N
                    ipt = jpt + itemp
                    if (ipt > N) then
                        itemp = jpt
                        jpt = ipt - N
                        ipt = itemp
                    end if
                    Xpt(Nf, ipt) = Xpt(ipt + 1, ipt)
                    Xpt(Nf, jpt) = Xpt(jpt + 1, jpt)
                elseif (nfm >= 1 .and. nfm <= N) then
                    stepa = Rhobeg
                    if (Su(nfm) == zero) stepa = -stepa
                    Xpt(Nf, nfm) = stepa
                elseif (nfm > N) then
                    stepa = Xpt(Nf - N, nfx)
                    stepb = -Rhobeg
                    if (Sl(nfx) == zero) stepb = DMIN1(two * Rhobeg, Su(&
     &nfx))
                    if (Su(nfx) == zero) stepb = DMAX1(-two * Rhobeg, Sl&
     &(nfx))
                    Xpt(Nf, nfx) = stepb
                end if
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
                do j = 1, N
                    X(j) = DMIN1(DMAX1(Xl(j), Xbase(j) + Xpt(Nf, j)), Xu&
     &(j))
                    if (Xpt(Nf, j) == Sl(j)) X(j) = Xl(j)
                    if (Xpt(Nf, j) == Su(j)) X(j) = Xu(j)
                end do
                call CALFUN(N, X, F)
                if (Iprint == 3) then
                    print 99001, Nf, F, (X(i), i=1, N)
      99001 format(/4X, 'Function number', I6, ' F =', 1PD18.10, ' The c&
     &orresponding X is:'/(2X, 5D15.6))
                end if
                Fval(Nf) = F
                if (Nf == 1) then
                    fbeg = F
                    Kopt = 1
                elseif (F < Fval(Kopt)) then
                    Kopt = Nf
                end if
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
                if (Nf > 2 * N + 1) then
                    ih = (ipt * (ipt - 1)) / 2 + jpt
                    Zmat(1, nfx) = recip
                    Zmat(Nf, nfx) = recip
                    Zmat(ipt + 1, nfx) = -recip
                    Zmat(jpt + 1, nfx) = -recip
                    temp = Xpt(Nf, ipt) * Xpt(Nf, jpt)
                    Hq(ih) = (fbeg - Fval(ipt + 1) - Fval(jpt + 1) + F) &
     &/ temp
                elseif (Nf >= 2 .and. Nf <= N + 1) then
                    Gopt(nfm) = (F - fbeg) / stepa
                    if (Npt < Nf + N) then
                        Bmat(1, nfm) = -one / stepa
                        Bmat(Nf, nfm) = one / stepa
                        Bmat(Npt + nfm, nfm) = -half * rhosq
                    end if
                elseif (Nf >= N + 2) then
                    ih = (nfx * (nfx + 1)) / 2
                    temp = (F - fbeg) / stepb
                    diff = stepb - stepa
                    Hq(ih) = two * (temp - Gopt(nfx)) / diff
                    Gopt(nfx) = (Gopt(nfx) * stepb - temp * stepa) / dif&
     &f
                    if (stepa * stepb < zero) then
                        if (F < Fval(Nf - N)) then
                            Fval(Nf) = Fval(Nf - N)
                            Fval(Nf - N) = F
                            if (Kopt == Nf) Kopt = Nf - N
                            Xpt(Nf - N, nfx) = stepb
                            Xpt(Nf, nfx) = stepa
                        end if
                    end if
                    Bmat(1, nfx) = -(stepa + stepb) / (stepa * stepb)
                    Bmat(Nf, nfx) = -half / Xpt(Nf - N, nfx)
                    Bmat(Nf - N, nfx) = -Bmat(1, nfx) - Bmat(Nf, nfx)
                    Zmat(1, nfx) = DSQRT(two) / (stepa * stepb)
                    Zmat(Nf, nfx) = DSQRT(half) / rhosq
                    Zmat(Nf - N, nfx) = -Zmat(1, nfx) - Zmat(Nf, nfx)
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
                end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     By Tom (on 04-06-2019):
!     If the evaluation returns an NaN or an infinity value, this
!     subroutine is stopped.
                if (F /= F .or. F > almost_infinity) exit
!     By Tom (on 04-06-2019):
!     If the target value is reached, stop the algorithm.
                if (F <= Ftarget) exit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if (Nf >= Npt .or. Nf >= Maxfun) exit
            end do
            end subroutine PRELIM