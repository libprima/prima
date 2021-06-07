!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the intersection-form version of lincoa.f90.
! The file is generated automatically and is NOT intended to be readable.
!
! In the intersection form, each continued line has an ampersand at column
! 73, and each continuation line has an ampersand at column 6. A Fortran
! file in such a form can be compiled both as fixed form and as free form.
!
! See http://fortranwiki.org/fortran/show/Continuation+lines for details.
!
! Generated using the interform.m script by Zaikun Zhang (www.zhangzk.net)
! on 07-Jun-2021.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*==lincoa.f90  processed by SPAG 7.50RE at 17:53 on 31 May 2021
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     1  MAXFUN,W)
            SUBROUTINE LINCOA(N,Npt,M,A,Ia,B,X,Rhobeg,Rhoend,Iprint,Maxf&
     &un,W, F,Info,Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IMPLICIT REAL*8*8 (A-H,O-Z)
            IMPLICIT NONE
!*--LINCOA12
!*++
!*++ Dummy argument declarations rewritten by SPAG
!*++
            INTEGER :: N
            INTEGER :: Npt
            INTEGER :: M
            REAL*8 , INTENT(IN) , DIMENSION(Ia,*) :: A
            INTEGER , INTENT(IN) :: Ia
            REAL*8 , INTENT(IN) , DIMENSION(*) :: B
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
            INTEGER :: i , iac , iamat , ib , ibmat , iflag , ifv , igo &
     &, ihq , ipq , ipqw , iqf , irc , irf , isp , istp , iw , ixb , ixn&
     & , ixo , ixp , ixs , izmat , j , ndim , np , nptm
            REAL*8 :: smallx , sum , temp , zero
!*++
!*++ End of declarations rewritten by SPAG
!*++
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     This subroutine seeks the least value of a function of many variables,
!       subject to general linear inequality constraints, by a trust region
!       method that forms quadratic models by interpolation. Usually there
!       is much freedom in each new model after satisfying the interpolation
!       conditions, which is taken up by minimizing the Frobenius norm of
!       the change to the second derivative matrix of the model. One new
!       function value is calculated on each iteration, usually at a point
!       where the current model predicts a reduction in the least value so
!       far of the objective function subject to the linear constraints.
!       Alternatively, a new vector of variables may be chosen to replace
!       an interpolation point that may be too far away for reliability, and
!       then the new point does not have to satisfy the linear constraints.
!       The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT must be set to the number of interpolation conditions, which is
!       required to be in the interval [N+2,(N+1)(N+2)/2]. Typical choices
!       of the author are NPT=N+6 and NPT=2*N+1. Larger values tend to be
!       highly inefficent when the number of variables is substantial, due
!       to the amount of work and extra difficulty of adjusting more points.
!     M must be set to the number of linear inequality constraints.
!     A is a matrix whose columns are the constraint gradients, which are
!       required to be nonzero.
!     IA is the first dimension of the array A, which must be at least N.
!     B is the vector of right hand sides of the constraints, the J-th
!       constraint being that the scalar product of A(.,J) with X(.) is at
!       most B(J). The initial vector X(.) is made feasible by increasing
!       the value of B(J) if necessary.
!     X is the vector of variables. Initial values of X(1),X(2),...,X(N)
!       must be supplied. If they do not satisfy the constraints, then B
!       is increased as mentioned above. X contains on return the variables
!       that have given the least calculated F subject to the constraints.
!     RHOBEG and RHOEND must be set to the initial and final values of a
!       trust region radius, so both must be positive with RHOEND<=RHOBEG.
!       Typically, RHOBEG should be about one tenth of the greatest expected
!       change to a variable, and RHOEND should indicate the accuracy that
!       is required in the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, the best
!       feasible vector of variables so far and the corresponding value of
!       the objective function are printed whenever RHO is reduced, where
!       RHO is the current lower bound on the trust region radius. Further,
!       each new value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN,
!       its value being at least NPT+1.
!     W is an array used for working space. Its length must be at least
!       M*(2+N) + NPT*(4+N+NPT) + N*(9+3*N) + MAX [ M+3*N, 2*M+N, 2*NPT ].
!       On return, W(1) is set to the final value of F, and W(2) is set to
!       the total number of function evaluations plus 0.5.
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
!     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
!       F to the value of the objective function for the variables X(1),
!       X(2),...,X(N). The value of the argument F is positive when CALFUN
!       is called if and only if the current X satisfies the constraints
!       to working accuracy.
!
!     Check that N, NPT and MAXFUN are acceptable.
!
            zero = 0.0D0
            smallx = 1.0D-6*Rhoend
            np = N + 1
            nptm = Npt - np
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      IF (N .LE. 1) THEN
!          IF (IPRINT .GT. 0) PRINT 10
!   10     FORMAT (/4X,'Return from LINCOA because N is less than 2.')
!          INFO=10
!          GOTO 80
!      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( Npt<N+2 .OR. Npt>((N+2)*np)/2 ) THEN
               IF ( Iprint>0 ) PRINT 99001
      99001 FORMAT (/4X,'Return from LINCOA because NPT is not in', ' th&
     &e required interval.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               Info = 5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               GOTO 99999
            ENDIF
            IF ( Maxfun<=Npt ) THEN
               IF ( Iprint>0 ) PRINT 99002
      99002 FORMAT (/4X,'Return from LINCOA because MAXFUN is less', ' t&
     &han NPT+1.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
               Info = 11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               GOTO 99999
            ENDIF
!
!     Normalize the constraints, and copy the resultant constraint matrix
!       and right hand sides into working space, after increasing the right
!       hand sides if necessary so that the starting point is feasible.
!
            iamat = MAX0(M+3*N,2*M+N,2*Npt) + 1
            ib = iamat + M*N
            iflag = 0
            IF ( M>0 ) THEN
               iw = iamat - 1
               DO j = 1 , M
                  sum = zero
                  temp = zero
                  DO i = 1 , N
                     sum = sum + A(i,j)*X(i)
                     temp = temp + A(i,j)**2
                  ENDDO
                  IF ( temp==zero ) THEN
                     IF ( Iprint>0 ) PRINT 99003
      99003 FORMAT (/4X,'Return from LINCOA because the gradient', ' of &
     &a constraint is zero.')
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                     Info = 12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     GOTO 99999
                  ENDIF
                  temp = DSQRT(temp)
                  IF ( sum-B(j)>smallx*temp ) iflag = 1
                  W(ib+j-1) = DMAX1(B(j),sum)/temp
                  DO i = 1 , N
                     iw = iw + 1
                     W(iw) = A(i,j)/temp
                  ENDDO
               ENDDO
            ENDIF
            IF ( iflag==1 ) THEN
               IF ( Iprint>0 ) PRINT 99004
      99004 FORMAT (/4X,'LINCOA has made the initial X feasible by', ' i&
     &ncreasing part(s) of B.')
            ENDIF
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
            ndim = Npt + N
            ixb = ib + M
            ixp = ixb + N
            ifv = ixp + N*Npt
            ixs = ifv + Npt
            ixo = ixs + N
            igo = ixo + N
            ihq = igo + N
            ipq = ihq + (N*np)/2
            ibmat = ipq + Npt
            izmat = ibmat + ndim*N
            istp = izmat + Npt*nptm
            isp = istp + N
            ixn = isp + Npt + Npt
            iac = ixn + N
            irc = iac + N
            iqf = irc + M
            irf = iqf + N*N
            ipqw = irf + (N*np)/2
!
!     The above settings provide a partition of W for subroutine LINCOB.
!
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
!     3  W(IRC),W(IQF),W(IRF),W(IPQW),W)
            CALL LINCOB(N,Npt,M,W(iamat),W(ib),X,Rhobeg,Rhoend,Iprint,Ma&
     &xfun, W(ixb),W(ixp),W(ifv),W(ixs),W(ixo),W(igo),W(ihq), W(ipq),W(i&
     &bmat),W(izmat),ndim,W(istp),W(isp),W(ixn), W(iac),W(irc),W(iqf),W(&
     &irf),W(ipqw),W,F,Info,Ftarget)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      99999 END SUBROUTINE LINCOA