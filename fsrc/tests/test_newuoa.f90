module test_solver_mod
! This module tests NEWUOA on a few simple problems.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Sunday, October 17, 2021 PM07:53:34

implicit none
private
public :: test_solver


contains


subroutine test_solver(probs, mindim, maxdim, dimstride, nrand)

use, non_intrinsic :: consts_mod, only : RP, IK, TWO, TEN, ZERO
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

character(len=PNLEN) :: probs_loc(100)
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
real(RP) :: rhoend
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: xhist(:, :)

if (present(probs)) then
    nprobs = int(size(probs), kind(nprobs))
    probs_loc(1:nprobs) = probs
else
    nprobs = 5_IK
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
        npt_list = [1_IK, &
            & n + 1_IK, n + 2_IK, n + 3_IK, &
            & 2_IK * n, 2_IK * n + 1_IK, 2_IK * n + 2_IK, &
            & (n + 1_IK) * (n + 2_IK) / 2_IK - 1_IK, (n + 1_IK) * (n + 2_IK) / 2_IK, &
            & (n + 1_IK) * (n + 2_IK) / 2_IK + 1_IK]
        do irand = 1, int(size(npt_list) + max(0_IK, nrand_loc), kind(irand))
            call setseed(int(sum(istr(probname)) + n + irand))  ! Initialize the random seed.
            if (irand <= size(npt_list)) then
                npt = npt_list(irand)
            else
                npt = int(floor(TEN * rand() * real(n, RP)), kind(npt))
            end if
            if (rand() <= 0.1_RP) then
                npt = 0
            end if
            maxfun = int(floor(2.0E2_RP * rand() * real(n, RP)), kind(maxfun))
            if (rand() <= 0.1_RP) then
                maxfun = 0
            end if
            maxhist = int(floor(TWO * rand() * real(maxfun, RP)), kind(maxhist))
            if (rand() <= 0.1_RP) then
                maxhist = 0
            end if
            rhobeg = getdelta0(n)
            rhoend = max(1.0E-6_RP, rhobeg * 1.0E1_RP**(6.0_RP * rand() - 5.0_RP))
            if (rand() <= 0.1_RP) then
                rhoend = rhobeg
            elseif (rand() <= 0.1_RP) then
                ! The probability to arrive here is 0.09. Note that the value of rand() changes.
                rhobeg = ZERO
            end if
            call safealloc(x, n) ! Not all compilers support automatic allocation yet, e.g., Absoft.
            x = getx0(n)
            print '(/1A, I3, 1A, I3)', trimstr(probname)//': N = ', n, ', Random test ', irand
            call newuoa(calfun, x, f, rhobeg=rhobeg, rhoend=rhoend, npt=npt, maxfun=maxfun, &
                & maxhist=maxhist, fhist=fhist, xhist=xhist, iprint=1_IK)
        end do
    end do
end do

end subroutine test_solver


end module test_solver_mod
