module test_solver_mod
!--------------------------------------------------------------------------------------------------!
! This module tests NEWUOA on a few simple problems.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Tuesday, December 21, 2021 AM12:12:52
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: test_solver


contains


subroutine test_solver(probs, mindim, maxdim, dimstride, nrand)

use, non_intrinsic :: consts_mod, only : RP, IK, TWO, TEN, ZERO, HUGENUM
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: cobyla_mod, only : cobyla
use, non_intrinsic :: noise_mod, only : noisy, noisy_calcfc, orig_calcfc
use, non_intrinsic :: param_mod, only : MINDIM_DFT, MAXDIM_DFT, DIMSTRIDE_DFT, NRAND_DFT
use, non_intrinsic :: prob_mod, only : PNLEN, problem_t, construct, destruct
use, non_intrinsic :: rand_mod, only : setseed, rand, randn
use, non_intrinsic :: string_mod, only : trimstr, istr

implicit none

character(len=PNLEN), intent(in), optional :: probs(:)
integer(IK), intent(in), optional :: mindim
integer(IK), intent(in), optional :: maxdim
integer(IK), intent(in), optional :: dimstride
integer(IK), intent(in), optional :: nrand

character(len=PNLEN) :: probname
character(len=PNLEN) :: probs_loc(100)
integer(IK) :: dimstride_loc
integer(IK) :: iprob
integer(IK) :: irand
integer(IK) :: m
integer(IK) :: maxdim_loc
integer(IK) :: maxfilt
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: mindim_loc
integer(IK) :: n
integer(IK) :: nprobs
integer(IK) :: nrand_loc
real(RP) :: cstrv
real(RP) :: ctol
real(RP) :: f
real(RP) :: ftarget
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP), allocatable :: chist(:)
real(RP), allocatable :: conhist(:, :)
real(RP), allocatable :: constr(:)
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)
type(problem_t) :: prob

if (present(probs)) then
    nprobs = int(size(probs), kind(nprobs))
    probs_loc(1:nprobs) = probs
else
    nprobs = 2_IK
    probs_loc(1:nprobs) = ['hexagon', 'hexagon']
end if

if (present(mindim)) then
    mindim_loc = mindim
else
    mindim_loc = MINDIM_DFT
end if

if (present(maxdim)) then
    maxdim_loc = maxdim
else
    maxdim_loc = MAXDIM_DFT
end if

if (present(dimstride)) then
    dimstride_loc = dimstride
else
    dimstride_loc = DIMSTRIDE_DFT
end if

if (present(nrand)) then
    nrand_loc = nrand
else
    nrand_loc = NRAND_DFT
end if

do iprob = 1, nprobs
    probname = probs_loc(iprob)
    call construct(prob, probname)  ! Construct the testing problem.
    m = prob % m
    n = prob % n
    do irand = 1, max(1_IK, nrand_loc)
        call setseed(int(sum(istr(probname)) + n + irand + IK + RP))  ! Initialize the random seed.
        maxfun = int(floor(2.0E2_RP * rand() * real(n, RP)), kind(maxfun))
        if (rand() <= 0.2_RP) then
            maxfun = 0
        end if
        maxhist = int(floor(TWO * rand() * real(maxfun, RP)), kind(maxhist))
        if (rand() <= 0.2_RP) then
            maxhist = 0
        end if
        maxfilt = int(floor(TWO * rand() * real(maxfun, RP)), kind(maxfilt))
        if (rand() <= 0.2_RP) then
            maxfilt = 0
        end if
        if (rand() <= 0.2_RP) then
            ctol = randn() * TEN**(-abs(TWO * randn()))
        elseif (rand() <= 0.2_RP) then  ! Note that the value of rand() changes.
            ctol = HUGENUM
        else
            ctol = ZERO
        end if
        if (rand() <= 0.2_RP) then
            ftarget = -TEN**abs(TWO * randn())
        elseif (rand() <= 0.2_RP) then  ! Note that the value of rand() changes.
            ftarget = HUGENUM
        else
            ftarget = -HUGENUM
        end if

        rhobeg = noisy(prob % Delta0)
        rhoend = max(1.0E-6_RP, rhobeg * 1.0E1_RP**(6.0_RP * rand() - 5.0_RP))
        if (rand() <= 0.2_RP) then
            rhoend = rhobeg
        elseif (rand() <= 0.2_RP) then  ! Note that the value of rand() changes.
            rhobeg = ZERO
        end if
        call safealloc(x, n) ! Not all compilers support automatic allocation yet, e.g., Absoft.
        x = noisy(prob % x0)
        orig_calcfc => prob % calcfc

        print '(/1A, I3, 1A, I3)', trimstr(probname)//': N = ', n, ', Random test ', irand
        call cobyla(noisy_calcfc, x, f, m=m, cstrv=cstrv, constr=constr, rhobeg=rhobeg, rhoend=rhoend, &
            & maxfun=maxfun, maxhist=maxhist, fhist=fhist, xhist=xhist, conhist=conhist, chist=chist, &
            & ctol=ctol, ftarget=ftarget, maxfilt=maxfilt, iprint=1_IK)

        deallocate (x)
        nullify (orig_calcfc)
    end do
    call destruct(prob)  ! Destruct the testing problem.
    ! DESTRUCT deallocates allocated arrays/pointers and nullify the pointers. Must be called.
end do

end subroutine test_solver


end module test_solver_mod
