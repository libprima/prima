! logging_mod is a module providing some subroutines for logging.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and papers.
!
! Last Modified: Wednesday, June 30, 2021 PM12:34:42


module logging_mod

implicit none
private
public :: logging


contains


subroutine logging(lnnmb, nf, f, nx, conv, mssg)
use consts_mod, only : RP, IK, MSSGLEN, OUTUNIT
use info_mod, only : FTARGET_ACHIEVED, MAXFUN_REACHED
use info_mod, only : SMALL_TR_RADIUS, TRSUBP_FAILED
use info_mod, only : NAN_X, NAN_INF_F, NAN_MODEL
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
character(len=100) :: fout
character(len=3) :: fstat
character(len=MSSGLEN) :: mssg
logical :: fexist


fout = solver//'_output.txt'
inquire (file=trim(fout), exist=fexist)
if (fexist) then
    fstat = 'old'
else
    fstat = 'new'
end if
open (unit=OUTUNIT, file=trim(fout), status=fstat, position='append', iostat=ios, action='write')
if (ios /= 0) then
    print '(1A)', 'Fail to open file '//trim(fout)//'!'
else
    write (OUTUNIT, '(/1A, I7, 4X, 1A, 1PD18.10, 4X, 1A, /(1P, 5D15.6))') 'Function number', nf, &
        & 'F = ', f, 'The corresponding X is:', x
    close (OUTUNIT)
end if

end subroutine fmssg


end module logging_mod
