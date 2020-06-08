! The gateway for UOBYQA
!
! Authors:
!     Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!     and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!     Department of Applied Mathematics,
!     The Hong Kong Polytechnic University.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

module fuobyqa
use pdfoconst ! See pdfoconst.F, which defines HUGENUM
implicit none
integer :: nf
double precision, allocatable :: fhist(:)
end module fuobyqa

subroutine muobyqa (n,x,rhobeg,rhoend,iprint,maxfun,w,f,info,funhist,ftarget)
use fuobyqa
implicit none
integer, intent(in) :: n,iprint,maxfun
integer, intent(out) :: info
double precision, intent(inout) :: x(n)
double precision, intent(in) :: rhobeg,rhoend,w((n*(42+n*(23+n*(8+n)))+max(2*n*n+4,18*n))/4+1),ftarget
double precision, intent(out) :: f,funhist(maxfun)

nf=0
if (allocated(fhist)) deallocate (fhist)
allocate(fhist(maxfun))
fhist(:)=hugenum

call uobyqa (n,x,rhobeg,rhoend,iprint,maxfun,w,f,info,ftarget)

funhist=fhist

deallocate(fhist)
return
end subroutine muobyqa

subroutine calfun (n,x,f)
use fuobyqa
implicit none
integer, intent(in) :: n
double precision, intent(in) :: x(n)
double precision, intent(out) :: f
double precision :: fun
external :: fun
f=fun(n,x)

nf=nf+1
fhist(nf)=f
return
end subroutine calfun
