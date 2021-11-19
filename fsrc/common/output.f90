module output_mod
!--------------------------------------------------------------------------------------------------!
! This module provides some subroutines concerning output to terminal/files. Note that these output
! operations are sequential in nature. In case parallelisum is desirable (especially during
! initializaton), the subroutines may have to be modified or disabled.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's Fortran 77 code and papers.
!
! Started: July 2020
!
! Last Modified: Friday, November 19, 2021 PM03:21:07
!--------------------------------------------------------------------------------------------------!

implicit none
private
public :: retmssg, rhomssg, fmssg

!--------------------------------------------------------------------------------------------------!
! Formats.
! Format for F at return: 16 digits for base, 4 digits for exponent.
character(len=*), parameter :: ffmt = '1PE25.16E4'
! Format for intermediate F: 10 digits for base, 4 digits for exponent.
character(len=*), parameter :: ffmt_intermediate = '1PE19.10E4'
! Format for X: 8 digits for base, 4 digits for exponent, 4 components in each line.
character(len=*), parameter :: xfmt = '/(1P, 4E19.8E4)'
! Format for RHO: 4 digits for base, 2 digits for exponent.
character(len=*), parameter :: rfmt = '1PE11.4E2'
! Format for integers: 10 digits.
character(len=*), parameter :: ifmt = 'I10'
! Separating spaces: 3 spaces.
character(len=*), parameter :: spaces = '3X'
! Format for NF at return.
character(len=*), parameter :: nf_fmt = '(1A, '//spaces//', 1A, '//ifmt//')'
! Format for F and X at return and when RHO is updated.
character(len=*), parameter :: fx_fmt = '(1A, '//ffmt//', '//spaces//', 1A, '//xfmt//')'
! Format for RHO and NF when RHO is updated.
character(len=*), parameter :: rnf_fmt = '(/1A, '//rfmt//', '//spaces//', 1A, '//ifmt//')'
! Format for NF, F, and X during iterations.
character(len=*), parameter :: nffx_fmt = &
    & '(/1A, '//ifmt//', '//spaces//', 1A, '//ffmt_intermediate//', '//spaces//', 1A, '//xfmt//')'
!--------------------------------------------------------------------------------------------------!


contains


subroutine retmssg(info, iprint, nf, f, x, solver)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages at return.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, MSSGLEN, OUTUNIT, FNAMELEN
use, non_intrinsic :: info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED
use, non_intrinsic :: info_mod, only : SMALL_TR_RADIUS, TRSUBP_FAILED
use, non_intrinsic :: info_mod, only : NAN_INF_X, NAN_INF_F, NAN_MODEL
implicit none

! Inputs
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: info
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)
character(len=*), intent(in) :: solver

! Local variables
integer :: ios  ! Should be an integer of default kind
character(len=FNAMELEN) :: fout
character(len=3) :: fstat
character(len=MSSGLEN) :: mssg
logical :: fexist


if (iprint == 0) then
    return
end if

if (info == FTARGET_ACHIEVED) then
    mssg = 'the target function value is achieved.'
else if (info == MAXFUN_REACHED) then
    mssg = 'the objective function has been evaluated MAXFUN times.'
else if (info == SMALL_TR_RADIUS) then
    mssg = 'the trust region radius reaches its lower bound.'
else if (info == TRSUBP_FAILED) then
    mssg = 'a trust region step has failed to reduce the quadratic model.'
else if (info == NAN_INF_X) then
    mssg = 'NaN or Inf occurs in x.'
else if (info == NAN_INF_F) then
    mssg = 'the objective function returns NaN/+Inf.'
else if (info == NAN_MODEL) then
    mssg = 'NaN occurs in the models.'
end if

if (iprint >= 1) then
    if (iprint >= 3) then
        print '(1X)'
    end if
    print '(/1A)', 'Return from '//solver//' because '//trim(mssg)
    print nf_fmt, 'At the return from '//solver, 'Number of function evaluations = ', nf
    print fx_fmt, 'Least function value = ', f, 'The corresponding X is:', x
    print '(1X)'
end if

if (iprint <= -1) then
    fout = solver//'_output.txt'
    inquire (file=trim(fout), exist=fexist)
    if (fexist) then
        fstat = 'old'
    else
        fstat = 'new'
    end if
    open (unit=OUTUNIT, file=trim(fout), status=fstat, position='append', iostat=ios, action='write')
    if (ios /= 0) then
        print '(/1A)', 'Fail to open file '//trim(fout)//'!'
    else
        if (iprint <= -3) then
            write (OUTUNIT, '(1X)')
        end if
        write (OUTUNIT, '(/1A)') 'Return from '//solver//' because '//trim(mssg)
        write (OUTUNIT, nf_fmt) 'At the return from '//solver, 'Number of function evaluations = ', nf
        write (OUTUNIT, fx_fmt) 'Least function value = ', f, 'The corresponding X is:', x
        close (OUTUNIT)
    end if
    !print '(/1A /)', 'The output is printed to ' // trim(fout) // '.'
end if

end subroutine retmssg


subroutine rhomssg(iprint, nf, f, rho, x, solver)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages when RHO is updated.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, OUTUNIT, FNAMELEN
implicit none

! Inputs
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: rho
real(RP), intent(in) :: x(:)
character(len=*), intent(in) :: solver

! Local variables
integer :: ios  ! Should be an integer of default kind
character(len=FNAMELEN) :: fout
character(len=3) :: fstat
logical :: fexist


if (iprint == 0) then
    return
end if

if (iprint >= 2) then
    if (iprint >= 3) then
        print '(1X)'
    end if
    print rnf_fmt, 'New RHO = ', rho, 'Number of function evaluations = ', nf
    print fx_fmt, 'Least function value = ', f, 'The corresponding X is:', x
end if

if (iprint <= -2) then
    fout = solver//'_output.txt'
    inquire (file=trim(fout), exist=fexist)
    if (fexist) then
        fstat = 'old'
    else
        fstat = 'new'
    end if
    open (unit=OUTUNIT, file=trim(fout), status=fstat, position='append', iostat=ios, action='write')
    if (ios /= 0) then
        print '(/1A)', 'Fail to open file '//trim(fout)//'!'
    else
        if (iprint <= -3) then
            write (OUTUNIT, '(1X)')
        end if
        write (OUTUNIT, rnf_fmt) 'New RHO = ', rho, 'Number of function evaluations = ', nf
        write (OUTUNIT, fx_fmt) 'Least function value = ', f, 'The corresponding X is:', x
        close (OUTUNIT)
    end if
end if

end subroutine rhomssg


subroutine fmssg(iprint, nf, f, x, solver)
!--------------------------------------------------------------------------------------------------!
! This subroutine prints messages at each iteration.
!--------------------------------------------------------------------------------------------------!
use, non_intrinsic :: consts_mod, only : RP, IK, OUTUNIT, FNAMELEN
implicit none

! Inputs
integer(IK), intent(in) :: iprint
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)
character(len=*), intent(in) :: solver

! Local variables
integer :: ios  ! Should be an integer of default kind
character(len=FNAMELEN) :: fout
character(len=3) :: fstat
logical :: fexist


if (iprint == 0) then
    return
end if

if (iprint >= 3) then
    print nffx_fmt, 'Function number', nf, 'F = ', f, 'The corresponding X is:', x
end if

if (iprint <= -3) then
    fout = solver//'_output.txt'
    inquire (file=trim(fout), exist=fexist)
    if (fexist) then
        fstat = 'old'
    else
        fstat = 'new'
    end if
    open (unit=OUTUNIT, file=trim(fout), status=fstat, position='append', iostat=ios, action='write')
    if (ios /= 0) then
        print '(/1A)', 'Fail to open file '//trim(fout)//'!'
    else
        write (OUTUNIT, nffx_fmt) 'Function number', nf, 'F = ', f, 'The corresponding X is:', x
        close (OUTUNIT)
    end if
end if

end subroutine fmssg


end module output_mod
