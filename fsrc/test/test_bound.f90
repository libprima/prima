program testbound

implicit none

integer :: iloop
integer :: m
integer :: mlist(1)
integer :: n
integer :: nloop

nloop = 2
do n = 1, 1
    do iloop = 1, nloop
        if (iloop <= size(mlist)) m = mlist(iloop)
    end do
end do

end program testbound
