subroutine biglag (n,npt,xopt,xpt,bmat,zmat,idz,ndim,knew, delta,d,alpha,hcol,gc,gd,s,w)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IMPLICIT REAL*8 (A-H,O-Z)
implicit real(kind(0.0d0)) (a-h,o-z)
implicit integer (i-n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dimension xopt(*),xpt(npt,*),bmat(ndim,*),zmat(npt,*),d(*), hcol(*),gc(*),gd(*),s(*),w(*)
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KNEW is the index of the interpolation point that is going to be moved.
!     DELTA is the current trust region bound.
!     D will be set to the step from XOPT to the new point.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     HCOL, GC, GD, S and W will be used for working space.
!
!     The step D is calculated in a way that attempts to maximize the modulus
!     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
!     the KNEW-th Lagrange function.
!
!     Set some constants.
!
half=0.5d0
one=1.0d0
zero=0.0d0
twopi=8.0d0*datan(one)
delsq=delta*delta
nptm=npt-n-1
!
!     Set the first NPT components of HCOL to the leading elements of the
!     KNEW-th column of H.
!
iterc=0
do k=1,npt
   hcol(k)=zero
end do
do j=1,nptm
   temp=zmat(knew,j)
   if (j < idz) temp=-temp
   do k=1,npt
hcol(k)=hcol(k)+temp*zmat(k,j)
   end do
end do
alpha=hcol(knew)
!
!     Set the unscaled initial direction D. Form the gradient of LFUNC at
!     XOPT, and multiply D by the second derivative matrix of LFUNC.
!
dd=zero
do i=1,n
   d(i)=xpt(knew,i)-xopt(i)
   gc(i)=bmat(knew,i)
   gd(i)=zero
   dd=dd+d(i)**2
end do
do k=1,npt
   temp=zero
   sum=zero
   do j=1,n
temp=temp+xpt(k,j)*xopt(j)
sum=sum+xpt(k,j)*d(j)
   end do
   temp=hcol(k)*temp
   sum=hcol(k)*sum
   do i=1,n
gc(i)=gc(i)+temp*xpt(k,i)
gd(i)=gd(i)+sum*xpt(k,i)
   end do
end do
!
!     Scale D and GD, with a sign change if required. Set S to another
!     vector in the initial two dimensional subspace.
!
gg=zero
sp=zero
dhd=zero
do i=1,n
   gg=gg+gc(i)**2
   sp=sp+d(i)*gc(i)
   dhd=dhd+d(i)*gd(i)
end do
scale=delta/dsqrt(dd)
if (sp*dhd < zero) scale=-scale
temp=zero
if (sp*sp > 0.99d0*dd*gg) temp=one
tau=scale*(dabs(sp)+half*scale*dabs(dhd))
if (gg*delsq < 0.01d0*tau*tau) temp=one
do i=1,n
   d(i)=scale*d(i)
   gd(i)=scale*gd(i)
   s(i)=gc(i)+temp*gd(i)
end do
!
!     Begin the iteration by overwriting S with a vector that has the
!     required length and direction, except that termination occurs if
!     the given D and S are nearly parallel.
!
   80 iterc=iterc+1
dd=zero
sp=zero
ss=zero
do i=1,n
   dd=dd+d(i)**2
   sp=sp+d(i)*s(i)
   ss=ss+s(i)**2
end do
temp=dd*ss-sp*sp
if (temp <= 1.0d-8*dd*ss) goto 160
denom=dsqrt(temp)
do i=1,n
   s(i)=(dd*s(i)-sp*d(i))/denom
   w(i)=zero
end do
!
!     Calculate the coefficients of the objective function on the circle,
!     beginning with the multiplication of S by the second derivative matrix.
!
do k=1,npt
   sum=zero
   do j=1,n
sum=sum+xpt(k,j)*s(j)
   end do
   sum=hcol(k)*sum
   do i=1,n
w(i)=w(i)+sum*xpt(k,i)
   end do
end do
cf1=zero
cf2=zero
cf3=zero
cf4=zero
cf5=zero
do i=1,n
   cf1=cf1+s(i)*w(i)
   cf2=cf2+d(i)*gc(i)
   cf3=cf3+s(i)*gc(i)
   cf4=cf4+d(i)*gd(i)
   cf5=cf5+s(i)*gd(i)
end do
cf1=half*cf1
cf4=half*cf4-cf1
!
!     Seek the value of the angle that maximizes the modulus of TAU.
!
taubeg=cf1+cf2+cf4
taumax=taubeg
tauold=taubeg
isave=0
iu=49
temp=twopi/dfloat(iu+1)
do i=1,iu
   angle=dfloat(i)*temp
   cth=dcos(angle)
   sth=dsin(angle)
   tau=cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth
   if (dabs(tau) > dabs(taumax)) then
taumax=tau
isave=i
tempa=tauold
   else if (i == isave+1) then
tempb=tau
   end if
   tauold=tau
end do
if (isave == 0) tempa=tau
if (isave == iu) tempb=taubeg
step=zero
if (tempa /= tempb) then
   tempa=tempa-taumax
   tempb=tempb-taumax
   step=half*(tempa-tempb)/(tempa+tempb)
end if
angle=temp*(dfloat(isave)+step)
!
!     Calculate the new D and GD. Then test for convergence.
!
cth=dcos(angle)
sth=dsin(angle)
tau=cf1+(cf2+cf4*cth)*cth+(cf3+cf5*cth)*sth
do i=1,n
   d(i)=cth*d(i)+sth*s(i)
   gd(i)=cth*gd(i)+sth*w(i)
   s(i)=gc(i)+gd(i)
end do
if (dabs(tau) <= 1.1d0*dabs(taubeg)) goto 160
if (iterc < n) goto 80
  160 return
end
