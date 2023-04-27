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
! Last Modified: Tuesday, May 03, 2022 PM10:02:07
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: update


contains


subroutine update(n, npt, xpt, bmat, zmat, idz, ndim, rsp, step, kopt, knew, vlag, w)

! Generic modules
use, non_intrinsic :: consts_mod, only : RP, IK, ONE, ZERO, HALF
use, non_intrinsic :: linalg_mod, only : matprod, inprod, planerot

implicit none

! Inputs
integer(IK), intent(in) :: kopt
integer(IK), intent(in) :: n
integer(IK), intent(in) :: ndim
integer(IK), intent(in) :: npt
real(RP), intent(in) :: rsp(2_IK * npt)
real(RP), intent(in) :: step(n)
real(RP), intent(in) :: xpt(n, npt)

! In-outputs
integer(IK), intent(inout) :: idz
integer(IK), intent(inout) :: knew
real(RP), intent(inout) :: bmat(n, npt + n)
real(RP), intent(inout) :: vlag(npt + n)
real(RP), intent(inout) :: w(npt + n)
real(RP), intent(inout) :: zmat(npt, npt - n - 1)

! Local variables
real(RP) :: alpha, beta, bsumm, denabs, denmax, denom, distsq,  &
&        dx, hdiag, scala, scalb, sqrtdn, ssq,  &
&        summ, tau, tausq, temp, tempa, tempb
integer(IK) :: i, iflag, j, ja, jb, jl, jp, k, nptm
real(RP) :: xopt(n)
real(RP) :: xdist(npt), vtmp(npt), xxpt(npt), sxpt(npt)!, wz(npt - n - 1), wzc(size(wz))
real(RP) :: grot(2, 2), ztest, zmatk1

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
nptm = npt - n - 1
!
!     Calculate VLAG and BETA for the current choice of STEP. The first NPT
!       elements of VLAG are set to the values of the Lagrange functions at
!       XPT(KOPT,.)+STEP(.). The first NPT compONEnts of W_check are held
!       in W, where W_check is defined in a paper on the updating method.
!

xxpt = matprod(xpt(:, kopt), xpt)
sxpt = matprod(step, xpt)
do k = 1, npt
    !w(k) = rsp(npt + k) * (HALF * rsp(npt + k) + rsp(k))
    w(k) = sxpt(k) * (HALF * sxpt(k) + xxpt(k))
    summ = ZERO
    do j = 1, n
        summ = summ + bmat(j, k) * step(j)
    end do
    !-----------------------------------------!
    !vlag(k) = summ
    vtmp(k) = summ
    !-----------------------------------------!
end do
beta = ZERO
!--------------------------------!
vlag(1:npt) = ZERO
!--------------------------------!
do k = 1, nptm
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
!------------------------------------------!
vlag(1:npt) = vlag(1:npt) + vtmp(1:npt)
!------------------------------------------!
bsumm = ZERO
dx = ZERO
ssq = ZERO
do j = 1, n
    summ = ZERO
    do i = 1, npt
        summ = summ + w(i) * bmat(j, i)
    end do
    bsumm = bsumm + summ * step(j)
    jp = npt + j
    do k = 1, n
        summ = summ + bmat(k, jp) * step(k)
    end do
    vlag(jp) = summ
    bsumm = bsumm + summ * step(j)
    dx = dx + step(j) * xpt(j, kopt)
    !ssq = ssq + step(j)**2
    ssq = ssq + step(j) * step(j)
end do

!!beta = dx * dx + ssq * (xxpt(kopt) + dx + dx + HALF * ssq) + beta - bsumm
!wz = matprod(w(1:npt), zmat); wzc = wz; wzc(1:idz - 1) = -wzc(1:idz - 1); beta = -inprod(wz, wzc)
!bsumm = sum(matprod(bmat(:, npt + 1:npt + n), step) * step(1:n) &
!     & + matprod(bmat(:, 1:npt), w(1:npt)) * step(1:n) &
!     & + matprod(bmat(:, 1:npt), w(1:npt)) * step(1:n))
!beta = dx**2 + ssq * (xxpt(kopt) + 2.0_RP * dx + HALF * ssq) + beta - bsumm
beta = inprod(vlag(1:npt + n), [w(1:npt), step])
xopt = xpt(:, kopt)
!beta = HALF * (inprod(xopt + step, xopt + step)**2 + inprod(xopt, xopt)**2) - inprod(xopt + step, xopt)**2 - beta
beta = dx**2 + ssq * (inprod(xopt, xopt) + dx + dx + half * ssq) &
    & - inprod(step, vlag(npt + 1:npt + n)) - inprod(w(1:npt), vlag(1:npt))
vlag(kopt) = vlag(kopt) + ONE

!
!     If KNEW is ZERO initially, then pick the index of the interpolation
!       point to be deleted, by maximizing the absolute value of the
!       denominator of the updating formula times a weighting factor.
!
!
if (knew == 0) then
    denmax = ZERO
    do k = 1, npt
        hdiag = ZERO
        do j = 1, nptm
            temp = ONE
            if (j < idz) temp = -ONE
            hdiag = hdiag + temp * zmat(k, j)**2
        end do
        hdiag = -sum(zmat(k, 1:idz - 1)**2) + sum(zmat(k, idz:size(zmat, 2))**2)
        denabs = abs(beta * hdiag + vlag(k)**2)
        distsq = ZERO
        do j = 1, n
            distsq = distsq + (xpt(j, k) - xpt(j, kopt))**2
        end do
        !temp = denabs * distsq * distsq
        temp = distsq**2 * denabs
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
    xdist = (sum((xpt - spread(xopt, dim=2, ncopies=npt))**2, dim=1))
    knew = maxloc(xdist, dim=1)
end if
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!
!     Apply the rotations that put ZEROs in the KNEW-th row of ZMAT.
!
ztest = 1.0E-20_RP * maxval(abs(zmat))
jl = 1
if (nptm >= 2) then
    do j = 2, nptm
        if (j == idz) then
            jl = idz
            !-----------------------------------------------------------------------------------------!
            ! Zaikun 20220410
            !else if (zmat(knew, j) /= ZERO) then
            !temp = sqrt(zmat(knew, jl)**2 + zmat(knew, j)**2)
            !tempa = zmat(knew, jl) / temp
            !tempb = zmat(knew, j) / temp
            !-------------------------------------------!
            ! Zaikun 20220412
            !else if (abs(zmat(knew, j)) > ZERO) then
        else if (abs(zmat(knew, j)) > ztest) then
            !-------------------------------------------!
            grot = planerot(zmat(knew, [jl, j]))  !!MATLAB: grot = planerot(zmat(knew, [jl, j])')
            tempa = grot(1, 1); tempb = grot(1, 2)
            !-----------------------------------------------------------------------------------------!
            do i = 1, npt
                temp = tempa * zmat(i, jl) + tempb * zmat(i, j)
                zmat(i, j) = tempa * zmat(i, j) - tempb * zmat(i, jl)
                zmat(i, jl) = temp
            end do
            zmat(knew, j) = ZERO
            !----------------------------------!
            !Zaikun 20220412
        else
            zmat(knew, j) = ZERO
            !----------------------------------!
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
tausq = tau**2
denom = alpha * beta + tausq
!--------------------------------------------------!
! Zaikun 20220410
! ?????????????????????????????
!if (denom == ZERO) then
!    knew = 0
!    goto 180
!end if
if (.not. abs(denom) > 0) return
!--------------------------------------------------!
sqrtdn = sqrt(abs(denom))
vlag(knew) = vlag(knew) - ONE
!
!     Complete the updating of ZMAT when there is only ONE nonZERO element
!       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to ONE when
!       the value of IDZ is going to be reduced.
!
iflag = 0
if (jl == 1) then
    tempa = tau / sqrtdn
    tempb = zmat(knew, 1) / sqrtdn
    zmatk1 = zmat(knew, 1)
    do i = 1, npt
        zmat(i, 1) = tempa * zmat(i, 1) - tempb * vlag(i)
        !zmat(i, 1) = (tau * zmat(i, 1) - zmatk1 * vlag(i)) / sqrtdn
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
    !scala = ONE / sqrt(abs(beta) * temp * temp + tausq)
    scala = ONE / sqrt(abs(beta) * temp**2 + tausq)
    scalb = scala * sqrtdn
    do i = 1, npt
        zmat(i, ja) = scala * (tau * zmat(i, ja) - temp * vlag(i))
        zmat(i, jb) = scalb * (zmat(i, jb) - tempa * w(i) - tempb * vlag(i))
    end do
    !----------------------------!
    !Zaikun 20220410
    !if (denom <= 0) then
    if (denom < 0) then
        !----------------------------!
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

function omega_inprod(idz, zmat, x, y) result(p)
!--------------------------------------------------------------------------------------------------!
! This function calculates P = X^T*OMEGA*Y, where, as Powell did in NEWUOA, BOBYQA, and LINCOA,
! OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
! OMEGA is the leading NPT-by-NPT block of the matrix H in (3.12) of the NEWUOA paper.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, DEBUGGING
use, non_intrinsic :: debug_mod, only : assert
use, non_intrinsic :: linalg_mod, only : matprod, inprod
implicit none

! Inputs
integer(IK), intent(in) :: idz
real(RP), intent(in) :: zmat(:, :)
real(RP), intent(in) :: x(:)
real(RP), intent(in) :: y(:)

! Outputs
real(RP) :: p

! Local variables
character(len=*), parameter :: srname = 'OMEGA_INPROD'
real(RP) :: xz(size(zmat, 2))
real(RP) :: yz(size(zmat, 2))

! Preconditions
if (DEBUGGING) then
    call assert(idz >= 1 .and. idz <= size(zmat, 2) + 1, '1 <= IDZ <= SIZE(ZMAT, 2) + 1', srname)
    call assert(size(x) == size(zmat, 1), 'SIZE(X) == SIZE(ZMAT, 1)', srname)
    call assert(size(y) == size(zmat, 1), 'SIZE(Y) == SIZE(ZMAT, 1)', srname)
end if

!====================!
! Calculation starts !
!====================!

xz = matprod(x, zmat)
xz(1:idz - 1) = -xz(1:idz - 1)
yz = matprod(y, zmat)
p = inprod(xz, yz)

!====================!
!  Calculation ends  !
!====================!

end function omega_inprod


end module update_mod
