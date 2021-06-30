!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of uobyqa.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 30-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==uobyqa.f90  processed by SPAG 7.50RE at 00:28 on 26 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      SUBROUTINE UOBYQA (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
            SUBROUTINE UOBYQA(N,X,Rhobeg,Rhoend,Iprint,Maxfun,W,F,Info, &
     &Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            IMPLICIT NONE
!*--UOBYQA11
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER :: N
            REAL*8 , DIMENSION(*) :: X
            REAL*8 :: Rhobeg
            REAL*8 , INTENT(INOUT) :: Rhoend
            INTEGER :: Iprint
            INTEGER :: Maxfun
            REAL*8 , DIMENSION(*) :: W
            REAL*8 :: F
            INTEGER :: Info
            REAL*8 :: Ftarget
!*++
!*++ Local variable declarations rewritten by SPAG
!*++
            INTEGER :: id , ig , ih , ipl , ipq , ivl , iw , ixb , ixn ,&
     & ixo , ixp , npt
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This subroutine seeks the least value of a function of many variables,
!     by a trust region method that forms quadratic models by interpolation.
!     The algorithm is described in "UOBYQA: unconstrained optimization by
!     quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA14,
!     University of Cambridge. The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
!       RHOBEG should be about one tenth of the greatest expected change to a
!       variable, and RHOEND should indicate the accuracy that is required in
!       the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       ( N**4 + 8*N**3 + 23*N**2 + 42*N + max [ 2*N**2 + 4, 18*N ] ) / 4.
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
!       -2: the objective function returns a NaN or nearly infinite
!           value.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
!     the value of the objective function for the variables X(1),X(2),...,X(N).
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
            npt = (N*N+3*N+2)/2
            ixb = 1
            ixo = ixb + N
            ixn = ixo + N
            ixp = ixn + N
            ipq = ixp + N*npt
            ipl = ipq + npt - 1
            ih = ipl + (npt-1)*npt
            ig = ih + N*N
            id = ig + N
            ivl = ih
            iw = id + N
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Zaikun, 2020-05-05
! When the data is passed from the interfaces to the Fortran code, RHOBEG,
! and RHOEND may change a bit (due to rounding ???). It was oberved in
! a MATLAB test that MEX passed 1 to Fortran as 0.99999999999999978.
! If we set RHOEND = RHOBEG in the interfaces, then it may happen
! that RHOEND > RHOBEG. That is why we do the following.
            Rhoend = MIN(Rhobeg,Rhoend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  W(IXN),W(IXP),W(IPQ),W(IPL),W(IH),W(IG),W(ID),W(IVL),W(IW))
            CALL UOBYQB(N,X,Rhobeg,Rhoend,Iprint,Maxfun,npt,W(ixb),W(ixo&
     &), W(ixn),W(ixp),W(ipq),W(ipl),W(ih),W(ig),W(id),W(ivl), W(iw),F,I&
     &nfo,Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            END SUBROUTINE UOBYQA