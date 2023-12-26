module recursive_mod
implicit none
private
public :: recursive_fun1

contains

subroutine chrosen(x, f, constr)
use, non_intrinsic :: consts_mod, only : RP
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
integer :: n
real(RP), parameter :: alpha = 4.0_RP
n = size(x)
f = sum((x(1:n - 1) - 1.0_RP)**2 + alpha * (x(2:n) - x(1:n - 1)**2)**2); 
constr = x(1:n - 1) - x(2:n)**2; 
end subroutine chrosen

subroutine recursive_fun1(x, f, constr)
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: cobyla_mod, only : cobyla
implicit none
real(RP), intent(in) :: x(:)
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
real(RP) :: x_loc(size(x))
x_loc = x
call cobyla(chrosen, int(size(x) - 1, IK), x_loc, f, nlconstr=constr)
end subroutine recursive_fun1

end module recursive_mod


module test_solver_mod
!--------------------------------------------------------------------------------------------------!
! This module tests COBYLA on a few simple problems.
!
! Coded by Zaikun ZHANG (www.zhangzk.net).
!
! Started: September 2021
!
! Last Modified: Wednesday, December 27, 2023 AM05:15:56
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: test_solver


contains


subroutine test_solver(probs, mindim, maxdim, dimstride, nrand, randseed, testdim)

use, non_intrinsic :: cobyla_mod, only : cobyla
use, non_intrinsic :: consts_mod, only : RP, IK, TWO, TEN, ZERO, REALMAX
use, non_intrinsic :: debug_mod, only : validate
use, non_intrinsic :: infnan_mod, only : is_neginf
use, non_intrinsic :: memory_mod, only : safealloc
use, non_intrinsic :: noise_mod, only : noisy, noisy_calcfc, orig_calcfc
use, non_intrinsic :: param_mod, only : MINDIM_DFT, MAXDIM_DFT, DIMSTRIDE_DFT, NRAND_DFT, RANDSEED_DFT
use, non_intrinsic :: prob_mod, only : PNLEN, PROB_T, construct, destruct
use, non_intrinsic :: rand_mod, only : setseed, rand, randn
use, non_intrinsic :: string_mod, only : strip, istr

implicit none

character(len=PNLEN), intent(in), optional :: probs(:)
integer(IK), intent(in), optional :: mindim
integer(IK), intent(in), optional :: maxdim
integer(IK), intent(in), optional :: dimstride
integer(IK), intent(in), optional :: nrand
integer, intent(in), optional :: randseed
character(len=*), intent(in), optional :: testdim

character(len=*), parameter :: bigprob = 'bigprob'
character(len=*), parameter :: solname = 'cobyla'
character(len=*), parameter :: srname = 'TEST_COBYLA'
character(len=:), allocatable :: testdim_loc
character(len=PNLEN) :: fix_dim_probs(100)  ! Problems with fixed dimensions
character(len=PNLEN) :: probname
character(len=PNLEN) :: probs_loc(100)  ! Maximal number of problems to test: 100
integer :: randseed_loc
integer :: rseed
integer(IK) :: dim_list(100)  ! Maximal number of dimensions to test: 100
integer(IK) :: dimstride_loc
integer(IK) :: idim
integer(IK) :: iprint
integer(IK) :: iprob
integer(IK) :: irand
integer(IK) :: m
integer(IK) :: maxdim_loc
integer(IK) :: maxfilt
integer(IK) :: maxfun
integer(IK) :: maxhist
integer(IK) :: mindim_loc
integer(IK) :: n
integer(IK) :: ndim
integer(IK) :: nprobs
integer(IK) :: nrand_loc
integer(IK), parameter :: bign = 80_IK
integer(IK), parameter :: largen = 1000_IK
real(RP) :: cstrv
real(RP) :: ctol
real(RP) :: f
real(RP) :: f_alt
real(RP) :: ftarget
real(RP) :: rhobeg
real(RP) :: rhoend
real(RP), allocatable :: Aineq(:, :)
real(RP), allocatable :: bineq(:)
real(RP), allocatable :: Aeq(:, :)
real(RP), allocatable :: beq(:)
real(RP), allocatable :: chist(:)
real(RP), allocatable :: nlchist(:, :)
real(RP), allocatable :: nlconstr(:)
real(RP), allocatable :: fhist(:)
real(RP), allocatable :: x(:)
real(RP), allocatable :: x0(:)
real(RP), allocatable :: x_alt(:)
real(RP), allocatable :: xhist(:, :)
real(RP), allocatable :: xl(:)
real(RP), allocatable :: xu(:)
type(PROB_T) :: prob

if (present(probs)) then
    nprobs = int(size(probs), kind(nprobs))
    probs_loc(1:nprobs) = probs
else
    nprobs = 12
    probs_loc(1:nprobs) = ['circle   ', 'ellipsoid', 'fletcheq1', 'fletcheq2', 'hs100    ', 'hexagon  ', 'rsnszk   ', &
        & 'chebyquad', 'chrosen  ', 'trigsabs ', 'trigssqs ', 'vardim   ']
end if
fix_dim_probs = '         '   ! Initialization, or compilers complain that the array is not (completely) defined.
fix_dim_probs(1:7) = ['circle   ', 'ellipsoid', 'fletcheq1', 'fletcheq2', 'hs100    ', 'hexagon  ', 'rsnszk   ']

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
    nrand_loc = NRAND_DFT * 5_IK  ! More random tests than default since we cannot vary NPT as other solvers.
end if

if (present(randseed)) then
    randseed_loc = randseed
else
    randseed_loc = RANDSEED_DFT
end if

if (present(testdim)) then
    testdim_loc = testdim
else
    testdim_loc = 'small'
end if


! Test the big problem
if (testdim_loc == 'big' .or. testdim_loc == 'large') then
    probname = bigprob
    n = merge(bign, largen, testdim_loc == 'big')
    do irand = 1, 1  ! The test is expensive
        rseed = int(sum(istr(solname)) + sum(istr(probname)) + n + irand + RP + randseed_loc)
        call setseed(rseed)
        m = int(min(int(10.0_RP * rand() * real(n, RP)), 10**min(range(0), range(0_IK))), IK)
        m = min(m, floor(real(huge(m)) / 8.0, IK) - n - 2_IK)  ! Avoid integer overflow when calculating UNIT_MEMO in PREPROC/HISTORY
        call construct(prob, probname, n, m)
        iprint = 2_IK
        if (int(n) + 2000 > huge(0_IK)) then
            maxfun = huge(0_IK)
        else
            maxfun = n + int(2000.0_RP * rand(), IK)
        end if
        maxhist = maxfun
        maxfilt = int(TWO * rand() * real(maxfun, RP), kind(maxfilt))
        if (rand() <= 0.5) then
            ctol = TEN**(-abs(5.0 * randn()))
        else
            ctol = ZERO
        end if
        ftarget = -REALMAX
        rhobeg = noisy(prob % Delta0)
        rhoend = max(1.0E-6_RP, rhobeg * 1.0E1_RP**(6.0_RP * rand() - 6.0_RP))
        call safealloc(x, n) ! Not all compilers support automatic allocation yet, e.g., Absoft.
        x = noisy(prob % x0)
        orig_calcfc => prob % calcfc
        call safealloc(xl, n)
        xl = prob % xl
        call safealloc(xu, n)
        xu = prob % xu
        call safealloc(Aineq, int(size(prob % Aineq, 1), IK), int(size(prob % Aineq, 2), IK))
        Aineq = prob % Aineq
        call safealloc(bineq, int(size(prob % bineq), IK))
        bineq = prob % bineq
        call safealloc(Aeq, int(size(prob % Aeq, 1), IK), int(size(prob % Aeq, 2), IK))
        Aeq = prob % Aeq
        call safealloc(beq, int(size(prob % beq), IK))
        beq = prob % beq

        print '(/A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)', &
           & strip(probname)//': N = ', n, ' M = ', m, ' Mineq = ', size(Aineq, 1), &
           & ' Meq = ', size(Aeq, 1), ', MAXFUN = ', maxfun, ', Random test ', irand

        call safealloc(nlconstr, m)
        call cobyla(noisy_calcfc, m, x, f, cstrv=cstrv, nlconstr=nlconstr, Aineq=Aineq, bineq=bineq, &
            & Aeq=Aeq, beq=beq, xl=xl, xu=xu, rhobeg=rhobeg, rhoend=rhoend, maxfun=maxfun, maxfilt=maxfilt,&
            & maxhist=maxhist, fhist=fhist, xhist=xhist, chist=chist, nlchist=nlchist,&
            & ftarget=ftarget, ctol=ctol, iprint=iprint)

        deallocate (x)
        deallocate (nlconstr)
        nullify (orig_calcfc)
        ! DESTRUCT deallocates allocated arrays/pointers and nullify the pointers. Must be called.
        call destruct(prob)  ! Destruct the testing problem.
    end do

else

    do iprob = 1, nprobs
        probname = probs_loc(iprob)
        if (any(probname == fix_dim_probs)) then
            call construct(prob, probname)  ! Construct the testing problem.
            ndim = 1
            dim_list(1) = prob % n
        else
            ndim = (maxdim_loc - mindim_loc) / dimstride_loc + 1_IK
            dim_list(1:ndim) = mindim_loc + dimstride_loc*[(idim - 1_IK, idim=1_IK, ndim)]
        end if
        do idim = 1, ndim
            if (any(probname == fix_dim_probs)) then
                call construct(prob, probname)
            else
                call construct(prob, probname, n=dim_list(idim))
            end if
            m = prob % m
            n = prob % n
            do irand = 1, max(1_IK, nrand_loc)
                ! Initialize the random seed using N, IRAND, RP, and RANDSEED_LOC. Do not include IK so
                ! that the results for different IK are the same.
                rseed = int(sum(istr(solname)) + sum(istr(probname)) + n + irand + RP + randseed_loc)
                call setseed(rseed)
                iprint = int(sign(min(3.0_RP, 1.5_RP * abs(randn())), randn()), kind(iprint))
                maxfun = int(2.0E2_RP * rand() * real(n, RP), kind(maxfun))
                if (rand() <= 0.2) then
                    maxfun = 0
                end if
                maxhist = int(TWO * rand() * real(max(10_IK * n, maxfun), RP), kind(maxhist))
                if (rand() <= 0.2) then
                    maxhist = -maxhist
                end if
                maxfilt = int(TWO * rand() * real(maxfun, RP), kind(maxfilt))
                if (rand() <= 0.2) then
                    maxfilt = 0
                end if
                if (rand() <= 0.2) then
                    ctol = randn() * TEN**(-abs(TWO * randn()))
                elseif (rand() <= 0.2) then  ! Note that the value of rand() changes.
                    ctol = REALMAX
                else
                    ctol = ZERO
                end if
                if (rand() <= 0.2) then
                    ftarget = -TEN**abs(TWO * randn())
                elseif (rand() <= 0.2) then  ! Note that the value of rand() changes.
                    ftarget = REALMAX
                else
                    ftarget = -REALMAX
                end if

                rhobeg = noisy(prob % Delta0)
                rhoend = max(1.0E-6_RP, rhobeg * 1.0E1_RP**(6.0_RP * rand() - 5.0_RP))
                if (rand() <= 0.2) then
                    rhoend = rhobeg
                elseif (rand() <= 0.2) then  ! Note that the value of rand() changes.
                    rhobeg = ZERO
                end if
                call safealloc(x0, n) ! Not all compilers support automatic allocation yet, e.g., Absoft.
                x0 = noisy(prob % x0)
                orig_calcfc => prob % calcfc
                call safealloc(xl, n)
                xl = prob % xl
                call safealloc(xu, n)
                xu = prob % xu
                call safealloc(Aineq, int(size(prob % Aineq, 1), IK), int(size(prob % Aineq, 2), IK))
                Aineq = prob % Aineq
                call safealloc(bineq, int(size(prob % bineq), IK))
                bineq = prob % bineq
                call safealloc(Aeq, int(size(prob % Aeq, 1), IK), int(size(prob % Aeq, 2), IK))
                Aeq = prob % Aeq
                call safealloc(beq, int(size(prob % beq), IK))
                beq = prob % beq

                print '(/A, I0, A, I0, A, I0, A, I0, A, I0, A, I0)', strip(probname)//': N = ', n, ' M = ', m, &
                    & ' Mineq = ', size(Aineq, 1), ' Meq = ', size(Aeq, 1), ', Random test ', irand

                call safealloc(x, n)
                x = x0
                call safealloc(nlconstr, m)
                call cobyla(noisy_calcfc, m, x, f, cstrv=cstrv, nlconstr=nlconstr, Aineq=Aineq, bineq=bineq, &
                    & Aeq=Aeq, beq=beq, xl=xl, xu=xu, rhobeg=rhobeg, rhoend=rhoend, &
                    & maxfun=maxfun, maxhist=maxhist, fhist=fhist, xhist=xhist, nlchist=nlchist, chist=chist, &
                    & ctol=ctol, ftarget=ftarget, maxfilt=maxfilt, iprint=iprint)

                if (prob % probtype == 'l') then  ! Run the test without nonlinear constraints
                    call safealloc(x_alt, n)
                    x_alt = x0
                    call cobyla(noisy_calcfc, m, x_alt, f_alt, Aineq=Aineq, bineq=bineq, Aeq=Aeq, beq=beq, &
                        & xl=xl, xu=xu, rhobeg=rhobeg, rhoend=rhoend, maxfun=maxfun, maxhist=maxhist, &
                        & fhist=fhist, xhist=xhist, ftarget=ftarget, maxfilt=maxfilt, iprint=iprint)
                    call validate(all(abs(x - x_alt) <= 0), 'X == X_ALT', srname)
                    call validate(abs(f - f_alt) <= 0 .or. (is_neginf(f) .and. is_neginf(f_alt)), 'F == F_ALT', srname)
                end if

                if (prob % probtype == 'b') then  ! Run the test without linear/nonlinear constraints
                    call safealloc(x_alt, n)
                    x_alt = x0
                    call cobyla(noisy_calcfc, m, x_alt, f_alt, xl=xl, xu=xu, rhobeg=rhobeg, rhoend=rhoend, &
                        & maxfun=maxfun, maxhist=maxhist, fhist=fhist, xhist=xhist, ftarget=ftarget, maxfilt=maxfilt, iprint=iprint)
                    call validate(all(abs(x - x_alt) <= 0), 'X == X_ALT', srname)
                    call validate(abs(f - f_alt) <= 0 .or. (is_neginf(f) .and. is_neginf(f_alt)), 'F == F_ALT', srname)
                end if

                if (prob % probtype == 'u') then  ! Run the test without constraints
                    call safealloc(x_alt, n)
                    x_alt = x0
                    call cobyla(noisy_calcfc, m, x_alt, f_alt, rhobeg=rhobeg, rhoend=rhoend, maxfun=maxfun, &
                        &maxhist=maxhist, fhist=fhist, xhist=xhist, ftarget=ftarget, maxfilt=maxfilt, iprint=iprint)
                    call validate(all(abs(x - x_alt) <= 0), 'X == X_ALT', srname)
                    call validate(abs(f - f_alt) <= 0 .or. (is_neginf(f) .and. is_neginf(f_alt)), 'F == F_ALT', srname)
                end if

                deallocate (x)
                deallocate (nlconstr)
                nullify (orig_calcfc)
            end do

            ! DESTRUCT deallocates allocated arrays/pointers and nullify the pointers. Must be called.
            call destruct(prob)  ! Destruct the testing problem.
        end do
    end do
end if

! Test recursive call.
! The depth of the recursion is 2. The first recursion is in RECURSIVE_FUN1, and the second is in
! RECURSIVE_FUN2. RECURSIVE_FUN1(Y) is defined by minimizing the CHROSEN function subject to
! X(1:N-1) - X(2:N)**2 <= 0 with Y being the starting point. RECURSIVE_FUN2(Y) is defined by
! RECURSIVE_FUN1 in a similar way. Note that RECURSIVE_FUN1 is essentially a constant function.
n = 3_IK
print '(/A, I0)', 'Testing recursive call of '//solname//' on a problem with N = ', n
call safealloc(x, n)
x = randn(n)
call safealloc(nlconstr, n - 1_IK)
call cobyla(recursive_fun2, n - 1_IK, x, f, nlconstr=nlconstr, iprint=2_IK)
deallocate (x, nlconstr)

contains

subroutine recursive_fun2(x_internal, f_internal, constr_internal)
use, non_intrinsic :: recursive_mod, only : recursive_fun1
implicit none
real(RP), intent(in) :: x_internal(:)
real(RP), intent(out) :: f_internal
real(RP), intent(out) :: constr_internal(:)
real(RP) :: x_loc(size(x_internal))
x_loc = x_internal
call cobyla(recursive_fun1, int(size(x_internal) - 1, IK), x_loc, f_internal, nlconstr=constr_internal)
end subroutine recursive_fun2

end subroutine test_solver


end module test_solver_mod
