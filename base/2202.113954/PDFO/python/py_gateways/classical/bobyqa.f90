! The gateway for BOBYQA
!
! Authors:
!     Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!     and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!     Department of Applied Mathematics,
!     The Hong Kong Polytechnic University.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

module fbobyqa
use pdfoconst ! See pdfoconst.F, which defines HUGENUM
implicit none
integer :: nf
double precision, allocatable :: fhist(:),chist(:),xlresmax(:),xuresmax(:)
end module fbobyqa

subroutine mbobyqa (n,npt,x,xl,xu,rhobeg,rhoend,iprint,maxfun,w,f,info,funhist,conhist,ftarget,resmax)
use fbobyqa
implicit none
integer, intent(in) :: n,npt,iprint,maxfun
integer, intent(out) :: info
integer :: i
double precision, intent(inout) :: x(n)
double precision, intent(in) :: xl(n),xu(n),rhobeg,rhoend,w((npt+5)*(npt+n)+3*n*(n+5)/2+1),ftarget
double precision, intent(out) :: f,funhist(maxfun),conhist(maxfun),resmax

nf=0
if (allocated(fhist)) deallocate (fhist)
allocate(fhist(maxfun))
fhist(:)=hugenum

if (allocated(chist)) deallocate (chist)
allocate(chist(maxfun))
chist(:)=hugenum

if (allocated(xlresmax)) deallocate (xlresmax)
allocate(xlresmax(n))
if (allocated(xuresmax)) deallocate (xuresmax)
allocate(xuresmax(n))
do i=1,n
    xlresmax(i)=xl(i)
    xuresmax(i)=xu(i)
enddo

call bobyqa (n,npt,x,xl,xu,rhobeg,rhoend,iprint,maxfun,w,f,info,ftarget)

funhist=fhist
conhist=chist

resmax=0.0d0
do i=1,n
    resmax=dmax1(resmax,x(i)-xu(i))
    resmax=dmax1(resmax,xl(i)-x(i))
enddo

deallocate(fhist)
deallocate(chist)
deallocate(xlresmax)
deallocate(xuresmax)
return
end subroutine mbobyqa

subroutine calfun (n,x,f)
use fbobyqa
implicit none
integer, intent(in) :: n
integer :: i
double precision, intent(in) :: x(n)
double precision, intent(out) :: f
double precision :: fun,resmax
external :: fun
f=fun(n,x)

resmax=0.0d0
do i=1,n
    resmax=dmax1(resmax,x(i)-xuresmax(i))
    resmax=dmax1(resmax,xlresmax(i)-x(i))
enddo

nf=nf+1
fhist(nf)=f
chist(nf)=resmax
return
end subroutine calfun
