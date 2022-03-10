! The gateway for NEWUOA
!
! Authors:
!     Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!     and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!     Department of Applied Mathematics,
!     The Hong Kong Polytechnic University.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

module fnewuoa
use pdfoconst ! See pdfoconst.F, which defines HUGENUM
implicit none
integer :: nf
double precision, allocatable :: fhist(:)
end module fnewuoa

subroutine mnewuoa (n,npt,x,rhobeg,rhoend,iprint,maxfun,w,f,info,funhist,ftarget)
use fnewuoa
implicit none
integer, intent(in) :: n,npt,iprint,maxfun
integer, intent(out) :: info
double precision, intent(inout) :: x(n)
double precision, intent(in) :: rhobeg,rhoend,w((npt+13)*(npt+n)+3*n*(n+3)/2+1),ftarget
double precision, intent(out) :: f,funhist(maxfun)

nf=0
if (allocated(fhist)) deallocate (fhist)
allocate(fhist(maxfun))
fhist(:)=hugenum

call newuoa (n,npt,x,rhobeg,rhoend,iprint,maxfun,w,f,info,ftarget)

funhist=fhist

deallocate(fhist)
return
end subroutine mnewuoa

subroutine calfun (n,x,f)
use fnewuoa
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
