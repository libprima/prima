module test_solver_mod
!--------------------------------------------------------------------------------------------------!
! This module tests NEWUOA on a few simple problems.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Wednesday, March 23, 2022 AM02:02:11
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: test_solver


contains


subroutine test_solver(probs, mindim, maxdim, dimstride, nrand, randseed)

use, non_intrinsic :: consts_mod, only : RP, IK, TWO, TEN, ZERO, HUGENUM
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: newuoa_mod, only : newuoa
use, non_intrinsic :: noise_mod, only : noisy, noisy_calfun, orig_calfun
use, non_intrinsic :: param_mod, only : MINDIM_DFT, MAXDIM_DFT, DIMSTRIDE_DFT, NRAND_DFT, RANDSEED_DFT
use, non_intrinsic :: prob_mod, only : PNLEN, PROB_T, construct, destruct
use, non_intrinsic :: rand_mod, only : setseed, rand, randn
use, non_intrinsic :: string_mod, only : trimstr, istr

implicit none

character(len=PNLEN), intent(in), optional :: probs(:)
integer(IK), intent(in), optional :: mindim
integer(IK), intent(in), optional :: maxdim
integer(IK), intent(in), optional :: dimstride
integer(IK), intent(in), optional :: nrand
integer, intent(in), optional :: randseed

character(len=*), parameter :: bigprob = 'bigprob'
character(len=PNLEN) :: probname
character(len=PNLEN) :: probs_loc(100)
integer :: randseed_loc
integer :: rseed
integer(IK), parameter :: bign = 1000_IK
integer(IK) :: dim_list(100)  ! Maximal number of dimensions to test: 100
integer(IK) :: dimstride_loc
integer(IK) :: idim
integer(IK) :: iprint
integer(IK) :: iprob
integer(IK) :: irand
integer(IK) :: maxdim_loc
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: mindim_loc
integer(IK) :: n
integer(IK) :: ndim
integer(IK) :: nprobs
integer(IK) :: npt
integer(IK) :: nnpt
integer(IK) :: npt_list(10)
integer(IK) :: nrand_loc
logical :: test_bigprob = .false.
real(RP) :: f
real(RP) :: ftarget
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)
type(PROB_T) :: prob

if (present(probs)) then
    nprobs = int(size(probs), kind(nprobs))
    probs_loc(1:nprobs) = probs
else
    nprobs = 5_IK
    probs_loc(1:nprobs) = ['chebyquad', 'chrosen  ', 'trigsabs ', 'trigssqs ', 'vardim   ']
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

if (present(randseed)) then
    randseed_loc = randseed
else
    randseed_loc = RANDSEED_DFT
end if


do iprob = 1, nprobs
    probname = probs_loc(iprob)
    ndim = (maxdim_loc - mindim_loc) / dimstride_loc + 1_IK
    dim_list(1:ndim) = mindim_loc + dimstride_loc*[(idim - 1_IK, idim=1_IK, ndim)]
    do idim = 1, ndim
        n = dim_list(idim)
        call construct(prob, probname, n)  ! Construct the testing problem.

        ! NPT_LIST defines some extreme values of NPT.
        nnpt = 10_IK
        npt_list(1:nnpt) = [1_IK, &
            & n + 1_IK, n + 2_IK, n + 3_IK, &
            & 2_IK * n, 2_IK * n + 1_IK, 2_IK * n + 2_IK, &
            & (n + 1_IK) * (n + 2_IK) / 2_IK - 1_IK, (n + 1_IK) * (n + 2_IK) / 2_IK, &
            & (n + 1_IK) * (n + 2_IK) / 2_IK + 1_IK]
        do irand = 1, nnpt + max(0_IK, nrand_loc)
            ! Initialize the random seed using N, IRAND, RP, and RANDSEED_LOC. Do not include IK so
            ! that the results for different IK are the same.
            rseed = int(sum(istr(probname)) + n + irand + RP + randseed_loc)
            call setseed(rseed)
            if (irand <= nnpt) then
                npt = npt_list(irand)
            else
                npt = int(TEN * rand() * real(n, RP), kind(npt))
            end if
            if (rand() <= 0.2_RP) then
                npt = 0
            end if
            iprint = int(sign(min(3.0_RP, 1.5_RP * abs(randn())), randn()), kind(iprint))
            maxfun = int(2.0E2_RP * rand() * real(n, RP), kind(maxfun))
            if (rand() <= 0.2_RP) then
                maxfun = 0
            end if
            maxhist = int(TWO * rand() * real(max(10_IK * n, maxfun), RP), kind(maxhist))
            if (rand() <= 0.2_RP) then
                maxhist = 0
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
            orig_calfun => prob % calfun

            print '(/1A, I0, 1A, I0)', trimstr(probname)//': N = ', n, ', Random test ', irand
            call newuoa(noisy_calfun, x, f, rhobeg=rhobeg, rhoend=rhoend, npt=npt, maxfun=maxfun, &
                & maxhist=maxhist, fhist=fhist, xhist=xhist, ftarget=ftarget, iprint=iprint)

            deallocate (x)
            nullify (orig_calfun)
        end do

        ! DESTRUCT deallocates allocated arrays/pointers and nullify the pointers. Must be called.
        call destruct(prob)  ! Destruct the testing problem.
    end do
end do


! Test the big problem
if (test_bigprob) then
    probname = bigprob
    n = bign
    call construct(prob, probname, n)
    nnpt = 2_IK
    npt_list(1:nnpt) = [2_IK * n + 1_IK, &
        & int(min(floor(real(10_IK**min(range(0), range(0_IK))) / 2.0), (int(n) + 1) * (int(n) + 2) / 2), IK)]
    do irand = 1, nnpt
        rseed = int(sum(istr(probname)) + n + irand + RP + randseed_loc)
        npt = npt_list(irand)
        iprint = 2_IK
        maxfun = int(minval([10**min(range(0), range(0_IK)), 10 * int(npt), int(npt) + 1000]), IK)
        maxhist = maxfun
        ftarget = -HUGENUM
        rhobeg = noisy(prob % Delta0)
        rhoend = max(1.0E-6_RP, rhobeg * 1.0E1_RP**(6.0_RP * rand() - 5.0_RP))
        call safealloc(x, n) ! Not all compilers support automatic allocation yet, e.g., Absoft.
        x = noisy(prob % x0)
        orig_calfun => prob % calfun

        print '(/1A, I0, 1A, I0)', trimstr(probname)//': N = ', n, ', Random test ', irand
        call newuoa(noisy_calfun, x, f, rhobeg=rhobeg, rhoend=rhoend, npt=npt, maxfun=maxfun, &
            & maxhist=maxhist, fhist=fhist, xhist=xhist, ftarget=ftarget, iprint=iprint)

        deallocate (x)
        nullify (orig_calfun)
    end do
    ! DESTRUCT deallocates allocated arrays/pointers and nullify the pointers. Must be called.
    call destruct(prob)  ! Destruct the testing problem.
end if


end subroutine test_solver


end module test_solver_mod
