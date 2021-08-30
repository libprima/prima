module evalfc_mod

implicit none
private
public :: evalfc


contains

subroutine evalfc(x, f, constr, cstrv)

! Generic modules
use consts_mod, only : RP, ZERO, HUGEFUN, HUGECON
use infnan_mod, only : is_nan

implicit none

! Input
real(RP), intent(in) :: x(:)

! Outputs
real(RP), intent(out) :: f
real(RP), intent(out) :: constr(:)
real(RP), intent(out) :: cstrv

if (any(is_nan(x))) then
    ! Set F and CONSTR to NaN. This is necessary if the initial X contains NaN.
    f = sum(x)
    constr = f
    cstrv = f
else
    call calcfc(size(x), size(constr), x, f, constr)  ! Evaluate F and CONSTR.

    ! Moderated extreme barrier: replace NaN/huge objective or constraint values with a large but
    ! finite value. This is naive, and better approaches surely exist.
    if (f > HUGEFUN .or. is_nan(f)) then
        f = HUGEFUN
    end if
    where (constr < -HUGECON .or. is_nan(constr))
        ! The constraint is CONSTR(X) >= 0, so NaN should be replaced with a large negative value.
        constr = -HUGECON  ! MATLAB code: constr(constr < -HUGECON | isnan(constr)) = -HUGECON
    end where

    ! Moderate huge positive values of CONSTR, or they may lead to Inf/NaN in subsequent calculations.
    ! This is NOT an extreme barrier.
    constr = min(HUGECON, constr)
    ! Moderate F similarly. (Shouldn't we?)
    f = max(-HUGEFUN, f)

    ! Evaluate the constraint violation for constraints CONSTR(X) >= 0.
    cstrv = maxval([-constr, ZERO])
end if

end subroutine evalfc

end module evalfc_mod
