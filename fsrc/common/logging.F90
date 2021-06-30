! logging_mod is a module providing some subroutines for logging.
!
! Coded by Zaikun Zhang in July 2020 based on Powell's Fortran 77 code and papers.
!
! Last Modified: Wednesday, June 30, 2021 PM02:45:38


module logging_mod

implicit none
private
public :: logging


contains


subroutine logging(srname, lnnum, nf, f, x, con, mssg)
use consts_mod, only : RP, IK, MSSGLEN
implicit none

! Inputs
integer(IK), optional, intent(in) :: lnnum
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), intent(in) :: x(:)
real(RP), optional, intent(in) :: con(:)
character(len=*), optional, intent(in) :: srname 
character(len=*), optional, intent(in) :: mssg

! Local variables
integer(IK), parameter :: LOGUNIT = 11
integer :: ios  ! Should be an integer of default kind
character(len=100) :: fout
character(len=3) :: fstat
logical :: fexist


fout = solver//'_log.txt'
inquire (file=trim(fout), exist=fexist)
if (fexist) then
    fstat = 'old'
else
    fstat = 'new'
end if
open (unit=LOGUNIT, file=trim(fout), status=fstat, position='append', iostat=ios, action='write')
if (ios /= 0) then
    print '(1A)', 'Fail to open file '//trim(fout)//'!'
else
    write (LOGUNIT, '(/1A, 1A, I7, 4X, 1A, 1PD18.10, 4X, 1A, /(1P, 5D15.6))') 'Line number', lnnum, 'Function number', nf, &
        & 'F = ', f, 'The corresponding X is:', x
    close (LOGUNIT)
end if

end subroutine fmssg


end module logging_mod
