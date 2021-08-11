!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of bobyqa.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 11-Aug-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==bobyqa.f90  processed by SPAG 7.50RE at 17:55 on 25 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  MAXFUN,W)
            subroutine BOBYQA(N, Npt, X, Xl, Xu, Rhobeg, Rhoend, Iprint,&
     & Maxfun, W, F, Info, Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            implicit none
!*--BOBYQA12
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
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
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            integer :: ibmat, id, ifv, igo, ihq, ipq, isl, isu, ivl, iw,&
     & ixa, ixb, ixn, ixo, ixp, izmat, j, jsl, jsu, ndim, np
            real*8 :: temp, zero
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This subroutine seeks the least value of a function of many variables,
!     by applying a trust region method that forms quadratic models by
!     interpolation. There is usually some freedom in the interpolation
!     conditions, which is taken up by minimizing the Frobenius norm of
!     the change to the second derivative of the model, beginning with the
!     zero matrix. The values of the variables are constrained by upper and
!     lower bounds. The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in
!       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
!       recommended.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
!       bounds, respectively, on X(I). The construction of quadratic models
!       requires XL(I) to be strictly less than XU(I) for each I. Further,
!       the contribution to a model from changes to the I-th variable is
!       damaged severely by rounding errors if XU(I)-XL(I) is too small.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND no greater than
!       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
!       expected change to a variable, while RHOEND should indicate the
!       accuracy that is required in the final values of the variables. An
!       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
!       is less than 2*RHOBEG.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     F is the objective function value when the algorithm exit.
!     INFO is the exit flag, which can be set to:
!       0: the lower bound for the trust region radius is reached.
!       1: the target function value is reached.
!       2: a trust region step has failed to reduce the quadratic model.
!       3: the objective function has been evaluated MAXFUN times.
!       4: much cancellation in a denominator.
!       5: NPT is not in the required interval.
!       6: one of the difference XU(I)-XL(I) is less than 2*RHOBEG.
!       7: rounding errors are becoming damaging.
!       8: rounding errors prevent reasonable changes to X.
!       9: the denominator of the updating formule is zero.
!       10: N should not be less than 2.
!       11: MAXFUN is less than NPT+1.
!       12: the gradient of constraint is zero.
!       -1: NaN occurs in x.
!       -2: the objective function returns a NaN or nearly infinite value.
!       -3: NaN occurs in BMAT or ZMAT.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
!     F to the value of the objective function for the current values of the
!     variables X(1),X(2),...,X(N), which are generated automatically in a
!     way that satisfies the bounds given in XL and XU.
!
!     Return if the value of NPT is unacceptable.
!
            np = N + 1
            if (Npt < N + 2 .or. Npt > ((N + 2) * np) / 2) then
                if (Iprint > 0) print 99001
99001 format(/4X, 'Return from BOBYQA because NPT is not in', ' the requ&
     &ired interval')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                Info = 5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                goto 99999
            end if
!
!     Partition the working space array, so that different parts of it can
!     be treated separately during the calculation of BOBYQB. The partition
!     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
!     space that is taken by the last array in the argument list of BOBYQB.
!
            ndim = Npt + N
            ixb = 1
            ixp = ixb + N
            ifv = ixp + N * Npt
            ixo = ifv + Npt
            igo = ixo + N
            ihq = igo + N
            ipq = ihq + (N * np) / 2
            ibmat = ipq + Npt
            izmat = ibmat + ndim * N
            isl = izmat + Npt * (Npt - np)
            isu = isl + N
            ixn = isu + N
            ixa = ixn + N
            id = ixa + N
            ivl = id + N
            iw = ivl + ndim
!
!     Return if there is insufficient space between the bounds. Modify the
!     initial X if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. The lower and upper
!     bounds on moves from the updated X are set now, in the ISL and ISU
!     partitions of W, in order to provide useful and exact information about
!     components of X that become within distance RHOBEG from their bounds.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun, 2020-05-05
! When the data is passed from the interfaces to the Fortran code, RHOBEG,
! XU and XL may change a bit (due to rounding ???). It was oberved in
! a MATLAB test that MEX passed 1 to Fortran as 0.99999999999999978.
! If we set RHOBEG = MIN(XU-XL)/2 in the interfaces, then it may happen
! that RHOBEG > MIN(XU-XL)/2. That is why we do the following. After
! this, INFO=6 should never occur.
            Rhobeg = min(0.5D0 * (1.0D0 - 1.0D-5) * minval(Xu(1:N) - Xl(&
     &1:N)), Rhobeg)
! For the same reason, we ensure RHOEND <= RHOBEG by the following.
            Rhoend = min(Rhobeg, Rhoend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            zero = 0.0D0
            do j = 1, N
                temp = Xu(j) - Xl(j)
                if (temp < Rhobeg + Rhobeg) then
                    if (Iprint > 0) print 99002
99002 format(/4X, 'Return from BOBYQA because one of the', ' differences&
     & XU(I)-XL(I)'/6X, ' is less than 2*RHOBEG.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                    Info = 6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    goto 99999
                end if
                jsl = isl + j - 1
                jsu = jsl + N
                W(jsl) = Xl(j) - X(j)
                W(jsu) = Xu(j) - X(j)
                if (W(jsl) >= -Rhobeg) then
                    if (W(jsl) >= zero) then
                        X(j) = Xl(j)
                        W(jsl) = zero
                        W(jsu) = temp
                    else
                        X(j) = Xl(j) + Rhobeg
                        W(jsl) = -Rhobeg
                        W(jsu) = DMAX1(Xu(j) - X(j), Rhobeg)
                    end if
                elseif (W(jsu) <= Rhobeg) then
                    if (W(jsu) <= zero) then
                        X(j) = Xu(j)
                        W(jsl) = -temp
                        W(jsu) = zero
                    else
                        X(j) = Xu(j) - Rhobeg
                        W(jsl) = DMIN1(Xl(j) - X(j), -Rhobeg)
                        W(jsu) = Rhobeg
                    end if
                end if
            end do
!
!     Make the call of BOBYQB.
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     2  NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW))
            call BOBYQB(N, Npt, X, Xl, Xu, Rhobeg, Rhoend, Iprint, Maxfu&
     &n, W(ixb), W(ixp), W(ifv), W(ixo), W(igo), W(ihq), W(ipq), W(ibmat&
     &), W(izmat), ndim, W(isl), W(isu), W(ixn), W(ixa), W(id), W(ivl) ,&
     & W(iw), F, Info, Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
99999 end subroutine BOBYQA