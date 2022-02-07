subroutine lincoa(n, npt, m, a, ia, b, x, rhobeg, rhoend, iprint, maxfun, w, f, info, ftarget)

implicit real(kind(0.0D0)) (a - h, o - z)
implicit integer(i - n)
dimension a(ia, *), b(*), x(*), w(*)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
smallx = 1.0D-6 * rhoend
np = n + 1
nptm = npt - np
if (npt < n + 2 .or. npt > ((n + 2) * np) / 2) then
    info = 5
    goto 80
end if
if (maxfun <= npt) then
    info = 11
    goto 80
end if
!
!     Normalize the constraints, and copy the resultant constraint matrix
!       and right hand sides into working space, after increasing the right
!       hand sides if necessary so that the starting point is feasible.
!
iamat = max0(m + 3 * n, 2 * m + n, 2 * npt) + 1
ib = iamat + m * n
iflag = 0
if (m > 0) then
    iw = iamat - 1
    do j = 1, m
        sum = zero
        temp = zero
        do i = 1, n
            sum = sum + a(i, j) * x(i)
            temp = temp + a(i, j)**2
        end do
        if (temp <= 0) then
            info = 12
            goto 80
        end if
        temp = dsqrt(temp)
        if (sum - b(j) > smallx * temp) iflag = 1
        w(ib + j - 1) = dmax1(b(j), sum) / temp
        do i = 1, n
            iw = iw + 1
            w(iw) = a(i, j) / temp
        end do
    end do
end if
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
ndim = npt + n
ixb = ib + m
ixp = ixb + n
ifv = ixp + n * npt
ixs = ifv + npt
ixo = ixs + n
igo = ixo + n
ihq = igo + n
ipq = ihq + (n * np) / 2
ibmat = ipq + npt
izmat = ibmat + ndim * n
istp = izmat + npt * nptm
isp = istp + n
ixn = isp + npt + npt
iac = ixn + n
irc = iac + n
iqf = irc + m
irf = iqf + n * n
ipqw = irf + (n * np) / 2
!
!     The above settings provide a partition of W for subroutine LINCOB.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun, 2020-05-05
! When the data is passed from the interfaces to the Fortran code, RHOBEG,
! and RHOEND may change a bit (due to rounding ???). It was oberved in
! a MATLAB test that MEX passed 1 to Fortran as 0.99999999999999978.
! If we set RHOEND = RHOBEG in the interfaces, then it may happen
! that RHOEND > RHOBEG. That is why we do the following.
rhoend = min(rhobeg, rhoend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call lincob(n, npt, m, w(iamat), w(ib), x, rhobeg, rhoend, iprint, &
& maxfun, w(ixb), w(ixp), w(ifv), w(ixs), w(ixo), w(igo), w(ihq), &
& w(ipq), w(ibmat), w(izmat), ndim, w(istp), w(isp), w(ixn), w(iac), &
& w(irc), w(iqf), w(irf), w(ipqw), w, f, info, ftarget)
80 return
end subroutine lincoa
