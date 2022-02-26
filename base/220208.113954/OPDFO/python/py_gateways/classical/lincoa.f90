! The gateway for LINCOA
!
! Authors:
!     Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!     and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!     Department of Applied Mathematics,
!     The Hong Kong Polytechnic University.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

module flincoa
use pdfoconst ! See pdfoconst.F, which defines HUGENUM
implicit none
integer :: nf,iaresmax,mresmax
double precision, allocatable :: fhist(:),chist(:),aresmax(:,:),bresmax(:)
end module flincoa

subroutine mlincoa (n,npt,m,a,ia,b,x,rhobeg,rhoend,iprint,maxfun,w,f,info,funhist,conhist,ftarget,resmax)
use flincoa
implicit none
integer, intent(in) :: n,npt,m,ia,iprint,maxfun
integer, intent(out) :: info
integer :: i,j
double precision, intent(inout) :: x(n)
double precision, intent(in) :: a(ia,m),b(m),rhobeg,rhoend,w(m*(2+n)+npt*(4+n+npt)+n*(9+3*n)+max(m+3*n,max(2*m+n,2*npt))),ftarget
double precision, intent(out) :: f,funhist(maxfun),conhist(maxfun),resmax
double precision :: cval

nf=0
if (allocated(fhist)) deallocate (fhist)
allocate(fhist(maxfun))
fhist(:)=hugenum

if (allocated(chist)) deallocate (chist)
allocate(chist(maxfun))
chist(:)=hugenum

if (allocated(aresmax)) deallocate (aresmax)
allocate(aresmax(ia,m))
if (allocated(bresmax)) deallocate (bresmax)
allocate(bresmax(m))
do j=1,m
    bresmax(j)=b(j)
    do i=1,ia
        aresmax(i,j)=a(i,j)
    enddo
enddo

iaresmax=ia
mresmax=m

call lincoa (n,npt,m,a,ia,b,x,rhobeg,rhoend,iprint,maxfun,w,f,info,ftarget)

funhist=fhist
conhist=chist

resmax=0.0d0
do j=1,m
    cval=-b(j)
    do i=1,ia
        cval=cval+a(i,j)*x(i)
    enddo
    if (cval .ne. cval) then
        resmax = cval ! Set resmax=NaN if constraint contains NaN
        exit
    else
        resmax=dmax1(resmax,cval)
    endif
enddo

deallocate(fhist)
deallocate(chist)
deallocate(aresmax)
deallocate(bresmax)
return
end subroutine mlincoa

subroutine calfun (n,x,f)
use flincoa
implicit none
integer, intent(in) :: n
integer :: i,j
double precision, intent(in) :: x(n)
double precision, intent(out) :: f
double precision :: fun,resmax,cval
external :: fun
f=fun(n,x)

resmax=0.0d0
do j=1,mresmax
    cval=-bresmax(j)
    do i=1,iaresmax
        cval=cval+aresmax(i,j)*x(i)
    enddo
    if (cval .ne. cval) then
        resmax = cval ! Set resmax=NaN if constraint contains NaN
        exit
    else
        resmax=dmax1(resmax,cval)
    endif
enddo

nf=nf+1
fhist(nf)=f
chist(nf)=resmax
return
end subroutine calfun
