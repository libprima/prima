module update_mod
!--------------------------------------------------------------------------------------------------!
! This module contains subroutines concerning the update of the interpolation set.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and the paper
!
! M. J. D. Powell, On fast trust region methods for quadratic models with linear constraints,
! Math. Program. Comput., 7:237--267, 2015
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
!
! Started: February 2022
!
! Last Modified: Wednesday, March 30, 2022 PM08:53:29
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: update


contains


subroutine update(kopt, rsp, step, xpt, idz, knew, bmat, zmat, vlag)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, HALF, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert, validate

implicit none

! Inputs
integer(IK), intent(in) :: kopt
real(RP), intent(in) :: rsp(:)  ! RSP(2*NPT)
real(RP), intent(in) :: step(:)  ! STEP(N)
real(RP), intent(in) :: xpt(:, :)  ! XPT(N, NPT)

! In-outputs
integer(IK), intent(inout) :: idz
integer(IK), intent(inout) :: knew
real(RP), intent(inout) :: bmat(:, :)  ! BMAT(N, NPT + N)
real(RP), intent(inout) :: zmat(:, :)  ! ZMAT(NPT, NPT - N - 1)

! Outputs
real(RP), intent(out) :: vlag(:)  ! VLAG(NPT + N)

! Local variables
character(len=*), parameter :: srname = 'UPDATE'
real(RP) :: w(size(vlag))
real(RP) :: alpha, beta, bsum, denabs, denmax, denom, distsq,  &
&        dx, hdiag, scala, scalb, sqrtdn, ssq,  &
&        summ, tau, tausq, temp, tempa, tempb
integer(IK) :: i, iflag, j, ja, jb, jl, jp, k
integer(IK) :: n
integer(IK) :: npt
real(RP) :: xopt(size(xpt, 1))
real(RP) :: xdist(size(xpt, 2))


! Sizes.
n = size(xpt, 1)
npt = size(xpt, 2)

! Preconditions
if (DEBUGGING) then
    call assert(n >= 1, 'N >= 1', srname)
    call assert(npt >= n + 2, 'NPT >= N+2', srname)
    call assert(knew >= 0 .and. knew <= npt, '0 <= KNEW <= NPT', srname)
    call assert(kopt >= 1 .and. kopt <= npt, '1 <= KOPT <= NPT', srname)
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(rsp) == 2_IK * npt, 'SIZE(RSP) == 2*NPT', srname)
    call assert(size(step) == n, 'SIZE(STEP) == N', srname)
    call assert(size(bmat, 1) == n .and. size(bmat, 2) == npt + n, 'SIZE(BMAT) == [N, NPT+N]', srname)
    call assert(size(zmat, 1) == npt .and. size(zmat, 2) == npt - n - 1_IK, 'SIZE(ZMAT) == [NPT, NPT-N-1]', srname)
    call assert(size(vlag) == npt + n, 'SIZE(VLAG) == NPT+N', srname)
end if

!
!     The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM ,SP and STEP are
!       identical to the corresponding arguments in SUBROUTINE LINCOB.
!     KOPT is such that XPT(KOPT,.) is the current trust region centre.
!     KNEW on exit is usually positive, and then it is the index of an
!       interpolation point to be moved to the position XPT(KOPT,.)+STEP(.).
!       It is set on entry either to its final value or to 0. In the latter
!       case, the final value of KNEW is chosen to maximize the denominator
!       of the matrix updating formula times a weighting factor.
!     VLAG and W are used for working space, the first NPT+N elements of
!       both of these vectors being required.
!
!     The arrays BMAT and ZMAT with IDZ are updated, the new matrices being
!       the ONEs that are suitable after the shift of the KNEW-th point to
!       the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to ZERO
!       occurs if the calculation fails due to a ZERO denominator in the
!       updating formula, which should never happen.
!
!     Set some constants.
!
!
!     Calculate VLAG and BETA for the current choice of STEP. The first NPT
!       elements of VLAG are set to the values of the Lagrange functions at
!       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
!       in W, where W_check is defined in a paper on the updating method.
!
do k = 1, npt
    w(k) = rsp(npt + k) * (HALF * rsp(npt + k) + rsp(k))
    summ = ZERO
    do j = 1, n
        summ = summ + bmat(j, k) * step(j)
    end do
    vlag(k) = summ
end do
beta = ZERO
do k = 1, npt - n - 1_IK
    summ = ZERO
    do i = 1, npt
        summ = summ + zmat(i, k) * w(i)
    end do
    if (k < idz) then
        beta = beta + summ * summ
        summ = -summ
    else
        beta = beta - summ * summ
    end if
    do i = 1, npt
        vlag(i) = vlag(i) + summ * zmat(i, k)
    end do
end do
bsum = ZERO
dx = ZERO
ssq = ZERO
do j = 1, n
    summ = ZERO
    do i = 1, npt
        summ = summ + w(i) * bmat(j, i)
    end do
    bsum = bsum + summ * step(j)
    jp = npt + j
    do k = 1, n
        summ = summ + bmat(k, jp) * step(k)
    end do
    vlag(jp) = summ
    bsum = bsum + summ * step(j)
    dx = dx + step(j) * xpt(j, kopt)
    ssq = ssq + step(j)**2
end do
beta = dx * dx + ssq * (rsp(kopt) + dx + dx + HALF * ssq) + beta - bsum
vlag(kopt) = vlag(kopt) + ONE
!
!     If KNEW is ZERO initially, then pick the index of the interpolation
!       point to be deleted, by maximizing the absolute value of the
!       denominator of the updating formula times a weighting factor.
!
!
if (knew == 0) then
    knew = 1  ! Without this, SIGSEV may occur due to uninitialized KNEW.
    denmax = ZERO
    do k = 1, npt
        hdiag = ZERO
        do j = 1, npt - n - 1_IK
            temp = ONE
            if (j < idz) temp = -ONE
            hdiag = hdiag + temp * zmat(k, j)**2
        end do
        denabs = abs(beta * hdiag + vlag(k)**2)
        distsq = ZERO
        do j = 1, n
            distsq = distsq + (xpt(j, k) - xpt(j, kopt))**2
        end do
        temp = denabs * distsq * distsq
        if (temp > denmax) then
            denmax = temp
            knew = k
        end if
    end do
end if

!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
! Zaikun 20220318: KNEW can be 0 due to NaN
if (knew == 0) then
    xopt = xpt(:, kopt)
    xdist = sqrt(sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1))
    ! MATLAB: xdist = sqrt(sum((xpt - xopt).^2))  % xopt should be a column!! Implicit expansion
    knew = maxloc(xdist, dim=1)
end if
call validate(1 <= knew .and. knew <= npt, '1 <= KNEW <= NPT', srname)
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!

!
!     Apply the rotations that put ZEROs in the KNEW-th row of ZMAT.
!
jl = 1
if (npt - n >= 3) then
    do j = 2, npt - n - 1_IK
        if (j == idz) then
            jl = idz
        else if (zmat(knew, j) /= ZERO) then
            temp = sqrt(zmat(knew, jl)**2 + zmat(knew, j)**2)
            ! Zaikun 20220304: TEMP can be 0 in single precision. Detected by Absoft and nvfortran. Should call GROT.
            tempa = zmat(knew, jl) / temp
            tempb = zmat(knew, j) / temp
            do i = 1, npt
                temp = tempa * zmat(i, jl) + tempb * zmat(i, j)
                zmat(i, j) = tempa * zmat(i, j) - tempb * zmat(i, jl)
                zmat(i, jl) = temp
            end do
            zmat(knew, j) = ZERO
        end if
    end do
end if
!
!     Put the first NPT compONEnts of the KNEW-th column of the Z Z^T matrix
!       into W, and calculate the parameters of the updating formula.
!
tempa = zmat(knew, 1)
if (idz >= 2) tempa = -tempa
if (jl > 1) tempb = zmat(knew, jl)
do i = 1, npt
    w(i) = tempa * zmat(i, 1)
    if (jl > 1) w(i) = w(i) + tempb * zmat(i, jl)
end do
alpha = w(knew)
tau = vlag(knew)
tausq = tau * tau
denom = alpha * beta + tausq
vlag(knew) = vlag(knew) - ONE
if (denom == ZERO) then
    knew = 0
    goto 180
end if
sqrtdn = sqrt(abs(denom))
!
!     Complete the updating of ZMAT when there is only ONE nonZERO element
!       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to ONE when
!       the value of IDZ is going to be reduced.
!
iflag = 0
if (jl == 1) then
    tempa = tau / sqrtdn
    tempb = zmat(knew, 1) / sqrtdn
    do i = 1, npt
        zmat(i, 1) = tempa * zmat(i, 1) - tempb * vlag(i)
    end do
    if (denom < 0) then
        if (idz == 1) then
            idz = 2
        else
            iflag = 1
        end if
    end if
else
!
!     Complete the updating of ZMAT in the alternative case.
!
    ja = 1
    if (beta >= 0) ja = jl
    jb = jl + 1 - ja
    temp = zmat(knew, jb) / denom
    tempa = temp * beta
    tempb = temp * tau
    temp = zmat(knew, ja)
    scala = ONE / sqrt(abs(beta) * temp * temp + tausq)
    scalb = scala * sqrtdn
    do i = 1, npt
        zmat(i, ja) = scala * (tau * zmat(i, ja) - temp * vlag(i))
        zmat(i, jb) = scalb * (zmat(i, jb) - tempa * w(i) - tempb * vlag(i))
    end do
    if (denom <= 0) then
        if (beta < 0) then
            idz = idz + 1
        else
            iflag = 1
        end if
    end if
end if
!
!     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
!       ZMAT^T factorization gains another positive element. Then exchange
!       the first and IDZ-th columns of ZMAT.
!
if (iflag == 1) then
! Zaikun 2020-06-28: I came here when reading the corresponding part of
! the NEWUOA code, which seems to have a bug. In NEWUOA, IDZ is redued
! only if IDZ >= 2, which is reasonable. Here there seems no such
! restriction. Why? Did Powell use a different definition for IDZ? In
! NEWUOA, IDZ is an intger used to represent the leading NPT sub-matrix
! of H, which is called OMEGA in the paper and represented in the code
! as
!
! OMEGA = sum_{K = 1}^{NPT-N-1} S_K*ZMAT(:, K)*ZMAT(:, K)',
! where S(1:IDZ-1) = -1 and S(IDZ:NPT-N-1) = 1.
!
! Indeed, theoretically, OMEGA is positive semidefinite, and S should be
! all positive. The negative entries of S result the rounding errors
! that cause OMEGA to lose the postive semidefiniteness. Therefore, in
! most cases, IDZ is small (e.g., IDZ=1, meaning that OMEGA has not lost
! the positive semidefiniteness), but it cannot be nonpositive in the
! NEWUOA code. Is it different in the LINCOA code??? To be studied.
! Unfortunately, Powell did not write a LINCOA paper!!!
!
! The BOBYQA code does not have this part --- it does not use IDZ at
! all. Why?
    idz = idz - 1
    do i = 1, npt
        temp = zmat(i, 1)
        zmat(i, 1) = zmat(i, idz)
        zmat(i, idz) = temp
    end do
end if
!
!     Finally, update the matrix BMAT.
!
do j = 1, n
    jp = npt + j
    w(jp) = bmat(j, knew)
    tempa = (alpha * vlag(jp) - tau * w(jp)) / denom
    tempb = (-beta * w(jp) - tau * vlag(jp)) / denom
    do i = 1, jp
        bmat(j, i) = bmat(j, i) + tempa * vlag(i) + tempb * w(i)
        if (i > npt) bmat(i - npt, jp) = bmat(j, i)
    end do
end do
180 return
end subroutine update


end module update_mod
