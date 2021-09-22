module test_newuoa_mod
! This module tests NEWUOA on a few simple problems.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Wednesday, September 22, 2021 PM11:48:31

implicit none
private
public :: test_newuoa


contains


subroutine test_newuoa(probs, mindim, maxdim, dimstride, nrand)

use, non_intrinsic :: consts_mod, only : RP, IK, TWO, TEN
use, non_intrinsic :: param_mod, only : MINDIM_DFT, MAXDIM_DFT, DIMSTRIDE_DFT, NRAND_DFT
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: newuoa_mod, only : newuoa
use, non_intrinsic :: prob_mod, only : PNLEN, probname, calfun, getx0, getdelta0
use, non_intrinsic :: rand_mod, only : setseed, rand
use, non_intrinsic :: string_mod, only : trimstr, istr

implicit none

character(len=PNLEN), optional, intent(in) :: probs(:)
integer(IK), optional, intent(in) :: mindim
integer(IK), optional, intent(in) :: maxdim
integer(IK), optional, intent(in) :: dimstride
integer(IK), optional, intent(in) :: nrand

character(len=PNLEN) :: probs_loc(1000)
integer(IK) :: dimstride_loc
integer(IK) :: iprob
integer(IK) :: irand
integer(IK) :: maxdim_loc
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: mindim_loc
integer(IK) :: n
integer(IK) :: nprobs
integer(IK) :: npt
integer(IK) :: npt_list(10)
integer(IK) :: nrand_loc
real(RP) :: f
real(RP) :: rhobeg
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)

if (present(probs)) then
    nprobs = int(size(probs), kind(nprobs))
    probs_loc(1:nprobs) = probs
else
    nprobs = 5
    probs_loc(1:nprobs) = ['chebyqad', 'chrosen ', 'trigsabs', 'trigssqs', 'vardim  ']
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
    do n = mindim_loc, maxdim_loc, dimstride_loc
        ! NPT_LIST defines some extreme values of NPT.
        npt_list = int([1, &
            & n + 1, n + 2, n + 3, &
            & 2 * n, 2 * n + 1, 2 * n + 2, &
            & (n + 1) * (n + 2) / 2 - 1, (n + 1) * (n + 2) / 2, (n + 1) * (n + 2) / 2 + 1], &
            & kind(npt))
        do irand = 1, int(size(npt_list) + max(0_IK, nrand_loc), kind(irand))
            call setseed(int(sum(istr(probname)) + n + irand))  ! Initialize the random seed.
            if (irand <= size(npt_list)) then
                npt = npt_list(irand)
            else
                npt = int(floor(TEN * rand() * real(n, RP)), kind(npt))
            end if
            maxfun = int(floor(5.0E2_RP * rand() * real(n, RP)), kind(maxfun))
            maxhist = int(floor(TWO * rand() * real(maxfun, RP)), kind(maxhist))
            rhobeg = getdelta0(n)
            call safealloc(x, n) ! Not all compilers support automatic allocation yet, e.g., Absoft.
            x = getx0(n)
            print '(/1A, I3, 1A, I3)', trimstr(probname)//': N = ', n, ', Random test ', irand
            call newuoa(calfun, x, f, rhobeg=rhobeg, npt=npt, maxfun=maxfun, maxhist=maxhist, &
                & fhist=fhist, xhist=xhist, iprint=1_IK)
        end do
    end do
end do

end subroutine test_newuoa


end module test_newuoa_mod
