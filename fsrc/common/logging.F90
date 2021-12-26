! logging_mod is a module providing some subroutines for logging.
!
! Coded by Zaikun ZHANG (www.zhangzk.net) in July 2020 based on Powell's Fortran 77 code and papers.
!
! Last Modified: Wednesday, September 01, 2021 AM10:58:33


module logging_mod

implicit none
private
public :: logging


contains


subroutine logging(logfile, srname, lnnum, nf, f, x, constr, conv, msg)
use, non_intrinsic :: consts_mod, only : RP, IK
implicit none

! Inputs
integer(IK), intent(in) :: lnnum
integer(IK), intent(in) :: nf
real(RP), intent(in) :: f
real(RP), optional, intent(in) :: x(:)
real(RP), optional, intent(in) :: constr(:)
real(RP), optional, intent(in) :: conv
character(len=*), intent(in) :: logfile
character(len=*), intent(in) :: srname
character(len=*), optional, intent(in) :: msg

! Local variables
integer, parameter :: LOGUNIT = 11
integer :: ios  ! Should be an integer of default kind
character(len=3) :: fstat
logical :: fexist


inquire (file=trim(logfile), exist=fexist)
if (fexist) then
    fstat = 'old'
else
    fstat = 'new'
end if
open (unit=LOGUNIT, file=trim(logfile), status=fstat, position='append', iostat=ios, action='write')
if (ios /= 0) then
    print '(1A)', 'Fail to open file '//trim(logfile)//'!'
else
    write (LOGUNIT, '(/1A, I5, 4X, 1A, I7, 4X, 1A, 1PD23.15)') 'In '//trim(srname)//', Line number', lnnum, &
        & 'Function number', nf, 'F = ', f
    if (present(constr)) then
        write (LOGUNIT, '(1A, 1PD23.15)') 'The constraint violation is:', conv
    end if
    if (present(msg)) then
        write (LOGUNIT, '(1A)') 'Message: '//msg
    end if
    if (present(x)) then
        write (LOGUNIT, '(1A, /(1P, 5D23.15))') 'The corresponding X is:', x
    end if
    if (present(constr)) then
        write (LOGUNIT, '(1A, /(1P, 5D23.15))') 'The constraint is:', constr
    end if
    close (LOGUNIT)
end if

end subroutine logging


end module logging_mod
