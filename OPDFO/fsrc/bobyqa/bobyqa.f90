subroutine bobyqa(n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, w, f, info, ftarget)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IMPLICIT REAL*8 (A-H,O-Z)
implicit real(kind(0.0D0)) (a - h, o - z)
implicit integer(i - n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dimension x(*), xl(*), xu(*), w(*)
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
np = n + 1
if (npt < n + 2 .or. npt > ((n + 2) * np) / 2) then
    info = 5
    go to 40
end if
!
!     Partition the working space array, so that different parts of it can
!     be treated separately during the calculation of BOBYQB. The partition
!     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
!     space that is taken by the last array in the argument list of BOBYQB.
!
ndim = npt + n
ixb = 1
ixp = ixb + n
ifv = ixp + n * npt
ixo = ifv + npt
igo = ixo + n
ihq = igo + n
ipq = ihq + (n * np) / 2
ibmat = ipq + npt
izmat = ibmat + ndim * n
isl = izmat + npt * (npt - np)
isu = isl + n
ixn = isu + n
ixa = ixn + n
id = ixa + n
ivl = id + n
iw = ivl + ndim
!
!     Return if there is insufficient space between the bounds. Modify the
!     initial X if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. The lower and upper
!     bounds on moves from the updated X are set now, in the ISL and ISU
!     partitions of W, in order to provide useful and exact information about
!     components of X that become within distance RHOBEG from their bounds.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun, 2020-05-05
! When the data is passed from the interfaces to the Fortran code, RHOBEG,
! XU and XL may change a bit (due to rounding ???). It was oberved in
! a MATLAB test that MEX passed 1 to Fortran as 0.99999999999999978.
! If we set RHOBEG = MIN(XU-XL)/2 in the interfaces, then it may happen
! that RHOBEG > MIN(XU-XL)/2. That is why we do the following. After
! this, INFO=6 should never occur.
rhobeg = min(0.5D0 * (1.0D0 - 1.0D-5) * minval(xu(1:n) - xl(1:n)), rhobeg)
! For the same reason, we ensure RHOEND <= RHOBEG by the following.
rhoend = min(rhobeg, rhoend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
zero = 0.0D0
do j = 1, n
    temp = xu(j) - xl(j)
    if (temp < rhobeg + rhobeg) then
        info = 6
        go to 40
    end if
    jsl = isl + j - 1
    jsu = jsl + n
    w(jsl) = xl(j) - x(j)
    w(jsu) = xu(j) - x(j)
    if (w(jsl) >= -rhobeg) then
        if (w(jsl) >= zero) then
            x(j) = xl(j)
            w(jsl) = zero
            w(jsu) = temp
        else
            x(j) = xl(j) + rhobeg
            w(jsl) = -rhobeg
            w(jsu) = dmax1(xu(j) - x(j), rhobeg)
        end if
    else if (w(jsu) <= rhobeg) then
        if (w(jsu) <= zero) then
            x(j) = xu(j)
            w(jsl) = -temp
            w(jsu) = zero
        else
            x(j) = xu(j) - rhobeg
            w(jsl) = dmin1(xl(j) - x(j), -rhobeg)
            w(jsu) = rhobeg
        end if
    end if
end do
!
!     Make the call of BOBYQB.
!
call bobyqb(n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, w(ixb), w(ixp), w(ifv), w(ixo), &
& w(igo), w(ihq), w(ipq), w(ibmat), w(izmat), &
& ndim, w(isl), w(isu), w(ixn), w(ixa), w(id), w(ivl), w(iw), f, info, &
& ftarget)
40 return
end
