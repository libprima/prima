module prob_mod
!--------------------------------------------------------------------------------------------------!
! This module implements the following testing problems.
!
! Unconstrained:
! chebyqad
! chrosen
! trigsabs
! trigssqs
! vardim
!
! Nonlinearly constrained:
! hexagon
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: pintrf_mod, only : FUN, FUNCON
implicit none

private
public :: PNLEN
public :: problem_t
public :: construct
public :: destruct

integer, parameter :: PNLEN = 64
integer, parameter :: SEED_DFT = 42  ! Default random seed used in trigsabs and trigssqs

type problem_t
    character(len=PNLEN) :: probname  ! Should be allocatable, which is not supported by Absoft 22.0
    character :: probtype
    integer(IK) :: m
    integer(IK) :: n
    real(RP), allocatable :: x0(:)
    real(RP), allocatable :: lb(:)
    real(RP), allocatable :: ub(:)
    real(RP), allocatable :: Aeq(:, :)
    real(RP), allocatable :: beq(:)
    real(RP), allocatable :: Aineq(:, :)
    real(RP), allocatable :: bineq(:)
    real(RP) :: Delta0
    procedure(FUN), nopass, pointer :: calfun => null()
    procedure(FUNCON), nopass, pointer :: calcfc => null()
end type problem_t


contains


subroutine construct(prob, probname, n)
!--------------------------------------------------------------------------------------------------!
! This subroutine constructs a derived type PROB of type PROBLEM_T.
! In F2003, this subroutine can be contained in the definition of PROBLEM_T, but the declaration of
! PROB must be changed to CLASS(PROBLEM_T).
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: string_mod, only : lower, trimstr
implicit none

! Inputs
character(len=*), intent(in) :: probname
integer(IK), intent(in), optional :: n

! Outputs
type(problem_t), intent(out) :: prob

! Local variables
character(len=*), parameter :: srname = 'CONSTRUCT'
integer(IK) :: n_loc

if (present(n)) then
    n_loc = n
else
    n_loc = 1
end if

select case (lower(trimstr(probname)))
case ('chebyqad')
    call construct_chebyqad(prob, n_loc)
case ('chrosen')
    call construct_chrosen(prob, n_loc)
case ('hexagon')
    call construct_hexagon(prob)
case ('trigsabs')
    call construct_trigsabs(prob, n_loc)
case ('trigssqs')
    call construct_trigssqs(prob, n_loc)
case ('vardim')
    call construct_vardim(prob, n_loc)
case default
    call errstop(srname, 'Unknown problem: '//trimstr(probname))
end select
end subroutine construct


subroutine destruct(prob)
!--------------------------------------------------------------------------------------------------!
! This subroutine destructs a derived type PROB of type PROBLEM_T.
! F2003 has the FINALIZATION mechanism for derived types, but not yet supported by all compilers.
! Here we code an explicit destructor. It must be called when PROB is not used anymore.
! In F2003, this subroutine can be contained in the definition of PROBLEM_T as a FINAL subroutine.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
type(problem_t), intent(inout) :: prob

!if (allocated(prob % probname)) then
!    deallocate (prob % probname)
!end if
if (allocated(prob % x0)) then
    deallocate (prob % x0)
end if
if (allocated(prob % lb)) then
    deallocate (prob % lb)
end if
if (allocated(prob % ub)) then
    deallocate (prob % ub)
end if
if (allocated(prob % Aeq)) then
    deallocate (prob % Aeq)
end if
if (allocated(prob % beq)) then
    deallocate (prob % beq)
end if
if (allocated(prob % Aineq)) then
    deallocate (prob % Aineq)
end if
if (allocated(prob % bineq)) then
    deallocate (prob % bineq)
end if
nullify (prob % calfun)
nullify (prob % calcfc)
end subroutine destruct


include 'chebyqad.f90'

include 'chrosen.f90'

include 'hexagon.f90'

include 'trigsabs.f90'

include 'trigssqs.f90'

include 'vardim.f90'

end module prob_mod
