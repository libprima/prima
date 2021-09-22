
program test_newuoa

implicit none

integer :: npt_list(10)
integer :: irand
integer :: npt

do irand = 1, 20
    npt_list = 1
    if (irand <= size(npt_list)) then
        npt = npt_list(irand)
    else
        npt = 0
    end if
    write (*, *) npt
end do

end program test_newuoa
