      subroutine trsapp(n, npt, tol, x, xpt, gq, hq, pq, delta, s,      &
     & crvmin, qred, info)
      
      use consts, only : rp, one, two, half, zero, pi
      use infnan
      use lina
      implicit none
      
      integer, intent(in) :: n, npt
      integer, intent(out) :: info
      real(kind = rp), intent(in) :: tol, delta, x(n), xpt(n, npt),     &
     & gq(n), hq((n*(n + 1))/2), pq(npt)
      real(kind = rp), intent(out) :: s(n), crvmin, qred

      real(kind = rp) :: d(n), g(n), hd(n), hs(n)

      call otrsapp (n,npt,x,transpose(xpt),gq,hq,pq,delta,s, d,g,hd,hs,crvmin)

      end subroutine trsapp
subroutine otrsapp (n,npt,xopt,xpt,gq,hq,pq,delta,step, d,g,hd,hs,crvmin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IMPLICIT REAL*8 (A-H,O-Z)
implicit real(kind(0.0d0)) (a-h,o-z)
implicit integer (i-n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dimension xopt(*),xpt(npt,*),gq(*),hq(*),pq(*),step(*), d(*),g(*),hd(*),hs(*)
!
!     N is the number of variables of a quadratic objective function, Q say.
!     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
!       in order to define the current quadratic model Q.
!     DELTA is the trust region radius, and has to be positive.
!     STEP will be set to the calculated trial step.
!     The arrays D, G, HD and HS will be used for working space.
!     CRVMIN will be set to the least curvature of H along the conjugate
!       directions that occur, except that it is set to zero if STEP goes
!       all the way to the trust region boundary.
!
!     The calculation of STEP begins with the truncated conjugate gradient
!     method. If the boundary of the trust region is reached, then further
!     changes to STEP may be made, each one being in the 2D space spanned
!     by the current STEP and the corresponding gradient of Q. Thus STEP
!     should provide a substantial reduction to Q within the trust region.
!
!     Initialization, which includes setting HD to H times XOPT.
!
half=0.5d0
zero=0.0d0
twopi=8.0d0*datan(1.0d0)
delsq=delta*delta
iterc=0
itermax=n
itersw=itermax
do i=1,n
   d(i)=xopt(i)
end do
goto 170
!
!     Prepare for the first line search.
!
   20 qred=zero
dd=zero
do i=1,n
   step(i)=zero
   hs(i)=zero
   g(i)=gq(i)+hd(i)
   d(i)=-g(i)
   dd=dd+d(i)**2
end do
crvmin=zero
if (dd == zero) goto 160
ds=zero
ss=zero
gg=dd
ggbeg=gg
!
!     Calculate the step to the trust region boundary and the product HD.
!
   40 iterc=iterc+1
temp=delsq-ss
bstep=temp/(ds+dsqrt(ds*ds+dd*temp))
goto 170
   50 dhd=zero
do j=1,n
   dhd=dhd+d(j)*hd(j)
end do
!
!     Update CRVMIN and set the step-length ALPHA.
!
alpha=bstep
if (dhd > zero) then
   temp=dhd/dd
   if (iterc == 1) crvmin=temp
   crvmin=dmin1(crvmin,temp)
   alpha=dmin1(alpha,gg/dhd)
end if
qadd=alpha*(gg-half*alpha*dhd)
qred=qred+qadd
!
!     Update STEP and HS.
!
ggsav=gg
gg=zero
do i=1,n
   step(i)=step(i)+alpha*d(i)
   hs(i)=hs(i)+alpha*hd(i)
   gg=gg+(g(i)+hs(i))**2
end do
!
!     Begin another conjugate direction iteration if required.
!
if (alpha < bstep) then
   if (qadd <= 0.01d0*qred) goto 160
   if (gg <= 1.0d-4*ggbeg) goto 160
   if (iterc == itermax) goto 160
   temp=gg/ggsav
   dd=zero
   ds=zero
   ss=zero
   do i=1,n
d(i)=temp*d(i)-g(i)-hs(i)
dd=dd+d(i)**2
ds=ds+d(i)*step(i)
ss=ss+step(i)**2
   end do
   if (ds <= zero) goto 160
   if (ss < delsq) goto 40
end if
crvmin=zero
itersw=iterc
!
!     Test whether an alternative iteration is required.
!
   90 if (gg <= 1.0d-4*ggbeg) goto 160
sg=zero
shs=zero
do i=1,n
   sg=sg+step(i)*g(i)
   shs=shs+step(i)*hs(i)
end do
sgk=sg+shs
angtest=sgk/dsqrt(gg*delsq)
if (angtest <= -0.99d0) goto 160
!
!     Begin the alternative iteration by calculating D and HD and some
!     scalar products.
!
iterc=iterc+1
temp=dsqrt(delsq*gg-sgk*sgk)
tempa=delsq/temp
tempb=sgk/temp
do i=1,n
   d(i)=tempa*(g(i)+hs(i))-tempb*step(i)
end do
goto 170
  120 dg=zero
dhd=zero
dhs=zero
do i=1,n
   dg=dg+d(i)*g(i)
   dhd=dhd+hd(i)*d(i)
   dhs=dhs+hd(i)*step(i)
end do
!
!     Seek the value of the angle that minimizes Q.
!
cf=half*(shs-dhd)
qbeg=sg+cf
qsav=qbeg
qmin=qbeg
isave=0
iu=49
temp=twopi/dfloat(iu+1)
do i=1,iu
   angle=dfloat(i)*temp
   cth=dcos(angle)
   sth=dsin(angle)
   qnew=(sg+cf*cth)*cth+(dg+dhs*cth)*sth
   if (qnew < qmin) then
qmin=qnew
isave=i
tempa=qsav
   else if (i == isave+1) then
tempb=qnew
   end if
   qsav=qnew
end do
if (isave == zero) tempa=qnew
if (isave == iu) tempb=qbeg
angle=zero
if (tempa /= tempb) then
   tempa=tempa-qmin
   tempb=tempb-qmin
   angle=half*(tempa-tempb)/(tempa+tempb)
end if
angle=temp*(dfloat(isave)+angle)
!
!     Calculate the new STEP and HS. Then test for convergence.
!
cth=dcos(angle)
sth=dsin(angle)
reduc=qbeg-(sg+cf*cth)*cth-(dg+dhs*cth)*sth
gg=zero
do i=1,n
   step(i)=cth*step(i)+sth*d(i)
   hs(i)=cth*hs(i)+sth*hd(i)
   gg=gg+(g(i)+hs(i))**2
end do
qred=qred+reduc
ratio=reduc/qred
if (iterc < itermax .and. ratio > 0.01d0) goto 90
  160 return
!
!     The following instructions act as a subroutine for setting the vector
!     HD to the vector D multiplied by the second derivative matrix of Q.
!     They are called from three different places, which are distinguished
!     by the value of ITERC.
!
  170 do i=1,n
   hd(i)=zero
end do
do k=1,npt
   temp=zero
   do j=1,n
temp=temp+xpt(k,j)*d(j)
   end do
   temp=temp*pq(k)
   do i=1,n
hd(i)=hd(i)+temp*xpt(k,i)
   end do
end do
ih=0
do j=1,n
   do i=1,j
ih=ih+1
if (i < j) hd(j)=hd(j)+hq(ih)*d(i)
hd(i)=hd(i)+hq(ih)*d(j)
   end do
end do
if (iterc == 0) goto 20
if (iterc <= itersw) goto 50
goto 120
end
