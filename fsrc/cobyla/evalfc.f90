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

    ! Handle "hidden constraints" (indicated by NaN/huge objective value or constraint violation)
    ! with a "moderated extreme barrier". Surely, this is naive, and better approaches exist.
    if (f > HUGEFUN .or. is_nan(f)) then
        f = HUGEFUN
    end if
    where (constr < -HUGECON .or. is_nan(constr))
        constr = -HUGECON  ! MATLAB code: constr(constr < -HUGECON | isnan(constr)) = -HUGECON
    end where

    ! Moderate excessively negative values of F and excessively positive values of CONSTR, or they
    ! may lead to Inf/NaN in subsequent calculations. This is NOT extreme barrier.
    f = max(-HUGEFUN, f)
    constr = min(HUGECON, constr)

    ! Evaluate the constraint violation for constraints CONSTR(X) >= 0.
    cstrv = maxval([-constr, ZERO])
end if

end subroutine evalfc

end module evalfc_mod
