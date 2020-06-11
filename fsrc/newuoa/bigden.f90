subroutine bigden (n,npt,xopt,xpt,bmat,zmat,idz,ndim,kopt, knew,d,w,vlag,beta,s,wvec,prod)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IMPLICIT REAL*8 (A-H,O-Z)
implicit real(kind(0.0d0)) (a-h,o-z)
implicit integer (i-n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dimension xopt(*),xpt(npt,*),bmat(ndim,*),zmat(npt,*),d(*), w(*),vlag(*),s(*),wvec(ndim,*),prod(ndim,*)
dimension den(9),denex(9),par(9)
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     D will be set to the step from XOPT to the new point, and on entry it
!       should be the D that was calculated by the last call of BIGLAG. The
!       length of the initial D provides a trust region bound on the final D.
!     W will be set to Wcheck for the final choice of D.
!     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
!     BETA will be set to the value that will occur in the updating formula
!       when the KNEW-th interpolation point is moved to its new position.
!     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
!       for working space.
!
!     D is calculated in a way that should provide a denominator with a large
!     modulus in the updating formula when the KNEW-th interpolation point is
!     shifted to the new position XOPT+D.
!
!     Set some constants.
!
half=0.5d0
one=1.0d0
quart=0.25d0
two=2.0d0
zero=0.0d0
twopi=8.0d0*datan(one)
nptm=npt-n-1
!
!     Store the first NPT elements of the KNEW-th column of H in W(N+1)
!     to W(N+NPT).
!
do k=1,npt
   w(n+k)=zero
end do
do j=1,nptm
   temp=zmat(knew,j)
   if (j < idz) temp=-temp
   do k=1,npt
w(n+k)=w(n+k)+temp*zmat(k,j)
   end do
end do
alpha=w(n+knew)
!
!     The initial search direction D is taken from the last call of BIGLAG,
!     and the initial S is set below, usually to the direction from X_OPT
!     to X_KNEW, but a different direction to an interpolation point may
!     be chosen, in order to prevent S from being nearly parallel to D.
!
dd=zero
ds=zero
ss=zero
xoptsq=zero
do i=1,n
   dd=dd+d(i)**2
   s(i)=xpt(knew,i)-xopt(i)
   ds=ds+d(i)*s(i)
   ss=ss+s(i)**2
   xoptsq=xoptsq+xopt(i)**2
end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
! Zaikun 2019-08-29: With the original code, if DS, DD, or SS is 
! NaN, KSAV will not get a value. This may cause Segmentation Fault
! because XPT(KSAV, :) will later be accessed. 
!      IF (DS*DS .GT. 0.99D0*DD*SS) THEN
if (.not. (ds*ds <= 0.99d0*dd*ss)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ksav=knew
   dtest=ds*ds/ss
   do k=1,npt
if (k /= kopt) then
    dstemp=zero
    sstemp=zero
    do i=1,n
  diff=xpt(k,i)-xopt(i)
  dstemp=dstemp+d(i)*diff
  sstemp=sstemp+diff*diff
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zaikun 2019-08-29: See the comments below line number 30
!              IF (DSTEMP*DSTEMP/SSTEMP .LT. DTEST) THEN
    if (.not. (dstemp*dstemp/sstemp >= dtest)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ksav=k
  dtest=dstemp*dstemp/sstemp
  ds=dstemp
  ss=sstemp
    end if
end if
   end do
   do i=1,n
s(i)=xpt(ksav,i)-xopt(i)
   end do
end if
ssden=dd*ss-ds*ds
iterc=0
densav=zero
!
!     Begin the iteration by overwriting S with a vector that has the
!     required length and direction.
!
   70 iterc=iterc+1
temp=one/dsqrt(ssden)
xoptd=zero
xopts=zero
do i=1,n
   s(i)=temp*(dd*s(i)-ds*d(i))
   xoptd=xoptd+xopt(i)*d(i)
   xopts=xopts+xopt(i)*s(i)
end do
!
!     Set the coefficients of the first two terms of BETA.
!
tempa=half*xoptd*xoptd
tempb=half*xopts*xopts
den(1)=dd*(xoptsq+half*dd)+tempa+tempb
den(2)=two*xoptd*dd
den(3)=two*xopts*dd
den(4)=tempa-tempb
den(5)=xoptd*xopts
do i=6,9
   den(i)=zero
end do
!
!     Put the coefficients of Wcheck in WVEC.
!
do k=1,npt
   tempa=zero
   tempb=zero
   tempc=zero
   do i=1,n
tempa=tempa+xpt(k,i)*d(i)
tempb=tempb+xpt(k,i)*s(i)
tempc=tempc+xpt(k,i)*xopt(i)
   end do
   wvec(k,1)=quart*(tempa*tempa+tempb*tempb)
   wvec(k,2)=tempa*tempc
   wvec(k,3)=tempb*tempc
   wvec(k,4)=quart*(tempa*tempa-tempb*tempb)
   wvec(k,5)=half*tempa*tempb
end do
do i=1,n
   ip=i+npt
   wvec(ip,1)=zero
   wvec(ip,2)=d(i)
   wvec(ip,3)=s(i)
   wvec(ip,4)=zero
   wvec(ip,5)=zero
end do
!
!     Put the coefficents of THETA*Wcheck in PROD.
!
do jc=1,5
   nw=npt
   if (jc == 2 .or. jc == 3) nw=ndim
   do k=1,npt
prod(k,jc)=zero
   end do
   do j=1,nptm
sum=zero
do k=1,npt
    sum=sum+zmat(k,j)*wvec(k,jc)
end do
if (j < idz) sum=-sum
do k=1,npt
    prod(k,jc)=prod(k,jc)+sum*zmat(k,j)
end do
   end do
   if (nw == ndim) then
do k=1,npt
    sum=zero
    do j=1,n
  sum=sum+bmat(k,j)*wvec(npt+j,jc)
    end do
    prod(k,jc)=prod(k,jc)+sum
end do
   end if
   do j=1,n
sum=zero
do i=1,nw
    sum=sum+bmat(i,j)*wvec(i,jc)
end do
prod(npt+j,jc)=sum
   end do
end do
!
!     Include in DEN the part of BETA that depends on THETA.
!
do k=1,ndim
   sum=zero
   do i=1,5
par(i)=half*prod(k,i)*wvec(k,i)
sum=sum+par(i)
   end do
   den(1)=den(1)-par(1)-sum
   tempa=prod(k,1)*wvec(k,2)+prod(k,2)*wvec(k,1)
   tempb=prod(k,2)*wvec(k,4)+prod(k,4)*wvec(k,2)
   tempc=prod(k,3)*wvec(k,5)+prod(k,5)*wvec(k,3)
   den(2)=den(2)-tempa-half*(tempb+tempc)
   den(6)=den(6)-half*(tempb-tempc)
   tempa=prod(k,1)*wvec(k,3)+prod(k,3)*wvec(k,1)
   tempb=prod(k,2)*wvec(k,5)+prod(k,5)*wvec(k,2)
   tempc=prod(k,3)*wvec(k,4)+prod(k,4)*wvec(k,3)
   den(3)=den(3)-tempa-half*(tempb-tempc)
   den(7)=den(7)-half*(tempb+tempc)
   tempa=prod(k,1)*wvec(k,4)+prod(k,4)*wvec(k,1)
   den(4)=den(4)-tempa-par(2)+par(3)
   tempa=prod(k,1)*wvec(k,5)+prod(k,5)*wvec(k,1)
   tempb=prod(k,2)*wvec(k,3)+prod(k,3)*wvec(k,2)
   den(5)=den(5)-tempa-half*tempb
   den(8)=den(8)-par(4)+par(5)
   tempa=prod(k,4)*wvec(k,5)+prod(k,5)*wvec(k,4)
   den(9)=den(9)-half*tempa
end do
!
!     Extend DEN so that it holds all the coefficients of DENOM.
!
sum=zero
do i=1,5
   par(i)=half*prod(knew,i)**2
   sum=sum+par(i)
end do
denex(1)=alpha*den(1)+par(1)+sum
tempa=two*prod(knew,1)*prod(knew,2)
tempb=prod(knew,2)*prod(knew,4)
tempc=prod(knew,3)*prod(knew,5)
denex(2)=alpha*den(2)+tempa+tempb+tempc
denex(6)=alpha*den(6)+tempb-tempc
tempa=two*prod(knew,1)*prod(knew,3)
tempb=prod(knew,2)*prod(knew,5)
tempc=prod(knew,3)*prod(knew,4)
denex(3)=alpha*den(3)+tempa+tempb-tempc
denex(7)=alpha*den(7)+tempb+tempc
tempa=two*prod(knew,1)*prod(knew,4)
denex(4)=alpha*den(4)+tempa+par(2)-par(3)
tempa=two*prod(knew,1)*prod(knew,5)
denex(5)=alpha*den(5)+tempa+prod(knew,2)*prod(knew,3)
denex(8)=alpha*den(8)+par(4)-par(5)
denex(9)=alpha*den(9)+prod(knew,4)*prod(knew,5)
!
!     Seek the value of the angle that maximizes the modulus of DENOM.
!
sum=denex(1)+denex(2)+denex(4)+denex(6)+denex(8)
denold=sum
denmax=sum
isave=0
iu=49
temp=twopi/dfloat(iu+1)
par(1)=one
do i=1,iu
   angle=dfloat(i)*temp
   par(2)=dcos(angle)
   par(3)=dsin(angle)
   do j=4,8,2
par(j)=par(2)*par(j-2)-par(3)*par(j-1)
par(j+1)=par(2)*par(j-1)+par(3)*par(j-2)
   end do
   sumold=sum
   sum=zero
   do j=1,9
sum=sum+denex(j)*par(j)
   end do
   if (dabs(sum) > dabs(denmax)) then
denmax=sum
isave=i
tempa=sumold
   else if (i == isave+1) then
tempb=sum
   end if
end do
if (isave == 0) tempa=sum
if (isave == iu) tempb=denold
step=zero
if (tempa /= tempb) then
   tempa=tempa-denmax
   tempb=tempb-denmax
   step=half*(tempa-tempb)/(tempa+tempb)
end if
angle=temp*(dfloat(isave)+step)
!
!     Calculate the new parameters of the denominator, the new VLAG vector
!     and the new D. Then test for convergence.
!
par(2)=dcos(angle)
par(3)=dsin(angle)
do j=4,8,2
   par(j)=par(2)*par(j-2)-par(3)*par(j-1)
   par(j+1)=par(2)*par(j-1)+par(3)*par(j-2)
end do
beta=zero
denmax=zero
do j=1,9
   beta=beta+den(j)*par(j)
   denmax=denmax+denex(j)*par(j)
end do
do k=1,ndim
   vlag(k)=zero
   do j=1,5
vlag(k)=vlag(k)+prod(k,j)*par(j)
   end do
end do
tau=vlag(knew)
dd=zero
tempa=zero
tempb=zero
do i=1,n
   d(i)=par(2)*d(i)+par(3)*s(i)
   w(i)=xopt(i)+d(i)
   dd=dd+d(i)**2
   tempa=tempa+d(i)*w(i)
   tempb=tempb+w(i)*w(i)
end do
if (iterc >= n) goto 340
if (iterc > 1) densav=dmax1(densav,denold)
if (dabs(denmax) <= 1.1d0*dabs(densav)) goto 340
densav=denmax
!
!     Set S to half the gradient of the denominator with respect to D.
!     Then branch for the next iteration.
!
do i=1,n
   temp=tempa*xopt(i)+tempb*d(i)-vlag(npt+i)
   s(i)=tau*bmat(knew,i)+alpha*temp
end do
do k=1,npt
   sum=zero
   do j=1,n
sum=sum+xpt(k,j)*w(j)
   end do
   temp=(tau*w(n+k)-alpha*vlag(k))*sum
   do i=1,n
s(i)=s(i)+temp*xpt(k,i)
   end do
end do
ss=zero
ds=zero
do i=1,n
   ss=ss+s(i)**2
   ds=ds+d(i)*s(i)
end do
ssden=dd*ss-ds*ds
if (ssden >= 1.0d-8*dd*ss) goto 70
!
!     Set the vector W before the RETURN from the subroutine.
!
  340 do k=1,ndim
   w(k)=zero
   do j=1,5
w(k)=w(k)+wvec(k,j)*par(j)
   end do
end do
vlag(kopt)=vlag(kopt)+one
return
end
