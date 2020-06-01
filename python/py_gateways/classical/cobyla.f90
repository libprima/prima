! The gateway for COBYLA
!
! Authors:
!     Tom M. RAGONNEAU (tom.ragonneau@connect.polyu.hk)
!     and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
!     Department of Applied Mathematics,
!     The Hong Kong Polytechnic University.
!
! Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

module fcobyla
use pdfoconst ! See pdfoconst.F, which defines HUGENUM
implicit none
integer :: nf
double precision, allocatable :: fhist(:),chist(:)
end module fcobyla

subroutine mcobyla (n,m,x,rhobeg,rhoend,iprint,maxfun,w,iact,f,info,funhist,conhist,ftarget,resmax,conval)
use fcobyla
implicit none
integer :: n,m,iprint,maxfun,iact(m+1),info
double precision :: x(n),rhobeg,rhoend,w(n*(3*n+2*m+11)+4*m+6),f,funhist(maxfun),conhist(maxfun),ftarget,resmax,conval(m)

nf=0
if (allocated(fhist)) deallocate (fhist)
allocate(fhist(maxfun))
if (allocated(chist)) deallocate (chist)
allocate(chist(maxfun))
fhist(:)=hugenum
chist(:)=hugenum

call cobyla (n,m,x,rhobeg,rhoend,iprint,maxfun,w,iact,f,info,ftarget,resmax,conval)

funhist=fhist
conhist=chist
deallocate(fhist)
deallocate(chist)
return
end subroutine mcobyla

subroutine calcfc (n,m,x,f,con)
use fcobyla
implicit none
integer :: n,m,i
double precision :: x(n),f,con(m),fun,resmax
external :: fun,confun
f=fun(n,x)

resmax=0.0d0
if (m .gt. 0) then
    ! The call to the constraint subroutine should be made only if a
    ! constraint function is supplied in the Python code. If m = 0,
    ! no such function is defined.
    call confun(n,m,x,con)
endif
do i=1,m
    resmax=dmax1(resmax,-con(i))
enddo

nf=nf+1
fhist(nf)=f
chist(nf)=resmax
return
end subroutine calcfc
