module prob_mod
!--------------------------------------------------------------------------------------------------!
! This module implements the following testing problems.
!
! Unconstrained:
! chebyquad
! chrosen
! trigsabs
! trigssqs
! vardim
!
! Linearly constrained:
! tetrahedron
!
! Nonlinearly constrained:
! circle
! ellipsoid
! fletcheq1
! fletcheq2
! hs100
! rsnszk
! hexagon
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK
use, non_intrinsic :: param_mod, only : RANDSEED_DFT
use, non_intrinsic :: pintrf_mod, only : OBJ, OBJCON
implicit none

private
public :: PNLEN
public :: PROB_T
public :: construct
public :: destruct

integer, parameter :: PNLEN = 64

type PROB_T
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
    procedure(OBJ), nopass, pointer :: calfun => null()
    procedure(OBJCON), nopass, pointer :: calcfc => null()
end type PROB_T


contains


subroutine construct(prob, probname, n)
!--------------------------------------------------------------------------------------------------!
! This subroutine constructs a derived type PROB of type PROB_T.
! In F2003, this subroutine can be contained in the definition of PROB_T, but the declaration of
! PROB must be changed to CLASS(PROB_T). See:
! https://fortran-lang.discourse.group/t/a-derived-type-containing-a-callback-function-as-a-member/2364/13
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : IK
use, non_intrinsic :: debug_mod, only : errstop
use, non_intrinsic :: string_mod, only : lower, trimstr
implicit none

! Inputs
character(len=*), intent(in) :: probname
integer(IK), intent(in), optional :: n

! Outputs
type(PROB_T), intent(out) :: prob

! Local variables
character(len=*), parameter :: srname = 'CONSTRUCT'
integer(IK) :: n_loc

if (present(n)) then
    n_loc = n
else
    n_loc = 1
end if

select case (lower(trimstr(probname)))
case ('chebyquad')
    call construct_chebyquad(prob, n_loc)
case ('chrosen')
    call construct_chrosen(prob, n_loc)
case ('circle')
    call construct_circle(prob)
case ('ellipsoid')
    call construct_ellipsoid(prob)
case ('fletcheq1')
    call construct_fletcheq1(prob)
case ('fletcheq2')
    call construct_fletcheq2(prob)
case ('hexagon')
    call construct_hexagon(prob)
case ('hs100')
    call construct_hs100(prob)
case ('rsnszk')
    call construct_rsnszk(prob)
case ('ptinsq')
    call construct_ptinsq(prob, n_loc)
case ('tetrahedron')
    call construct_tetrahedron(prob)
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
! This subroutine destructs a derived type PROB of type PROB_T.
! F2003 has the FINALIZATION mechanism for derived types, but not yet supported by all compilers.
! Here we code an explicit destructor. It must be called when PROB is not used anymore.
! In F2003, this subroutine can be contained in the definition of PROB_T as a FINAL subroutine.
!--------------------------------------------------------------------------------------------------!
implicit none
! Inputs
type(PROB_T), intent(inout) :: prob

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


include 'chebyquad.f90'

include 'chrosen.f90'

include 'circle.f90'

include 'ellipsoid.f90'

include 'fletcheq1.f90'

include 'fletcheq2.f90'

include 'hexagon.f90'

include 'hs100.f90'

include 'rsnszk.f90'

include 'ptinsq.f90'

include 'tetrahedron.f90'

include 'trigsabs.f90'

include 'trigssqs.f90'

include 'vardim.f90'

end module prob_mod
