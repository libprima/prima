subroutine update (n,npt,bmat,zmat,idz,ndim,vlag,beta,knew,w)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IMPLICIT REAL*8 (A-H,O-Z)
implicit real(kind(0.0d0)) (a-h,o-z)
implicit integer (i-n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dimension bmat(ndim,*),zmat(npt,*),vlag(*),w(*)
!
!     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
!     interpolation point that has index KNEW. On entry, VLAG contains the
!     components of the vector Theta*Wcheck+e_b of the updating formula
!     (6.11), and BETA holds the value of the parameter that has this name.
!     The vector W is used for working space.
!
!     Set some constants.
!
one=1.0d0
zero=0.0d0
nptm=npt-n-1
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
jl=1
do j=2,nptm
   if (j == idz) then
jl=idz
   else if (zmat(knew,j) /= zero) then
temp=dsqrt(zmat(knew,jl)**2+zmat(knew,j)**2)
tempa=zmat(knew,jl)/temp
tempb=zmat(knew,j)/temp
do i=1,npt
    temp=tempa*zmat(i,jl)+tempb*zmat(i,j)
    zmat(i,j)=tempa*zmat(i,j)-tempb*zmat(i,jl)
    zmat(i,jl)=temp
end do
zmat(knew,j)=zero
   end if
end do
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
tempa=zmat(knew,1)
if (idz >= 2) tempa=-tempa
if (jl > 1) tempb=zmat(knew,jl)
do i=1,npt
   w(i)=tempa*zmat(i,1)
   if (jl > 1) w(i)=w(i)+tempb*zmat(i,jl)
end do
alpha=w(knew)
tau=vlag(knew)
tausq=tau*tau
denom=alpha*beta+tausq
vlag(knew)=vlag(knew)-one
!
!     Complete the updating of ZMAT when there is only one nonzero element
!     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
!     then the first column of ZMAT will be exchanged with another one later.
!
iflag=0
if (jl == 1) then
   temp=dsqrt(dabs(denom))
   tempb=tempa/temp
   tempa=tau/temp
   do i=1,npt
zmat(i,1)=tempa*zmat(i,1)-tempb*vlag(i)
   end do
   if (idz == 1 .and. temp < zero) idz=2
   if (idz >= 2 .and. temp >= zero) iflag=1
else
!
!     Complete the updating of ZMAT in the alternative case.
!
   ja=1
   if (beta >= zero) ja=jl
   jb=jl+1-ja
   temp=zmat(knew,jb)/denom
   tempa=temp*beta
   tempb=temp*tau
   temp=zmat(knew,ja)
   scala=one/dsqrt(dabs(beta)*temp*temp+tausq)
   scalb=scala*dsqrt(dabs(denom))
   do i=1,npt
zmat(i,ja)=scala*(tau*zmat(i,ja)-temp*vlag(i))
zmat(i,jb)=scalb*(zmat(i,jb)-tempa*w(i)-tempb*vlag(i))
   end do
   if (denom <= zero) then
if (beta < zero) idz=idz+1
if (beta >= zero) iflag=1
   end if
end if
!
!     IDZ is reduced in the following case, and usually the first column
!     of ZMAT is exchanged with a later one.
!
if (iflag == 1) then
   idz=idz-1
   do i=1,npt
temp=zmat(i,1)
zmat(i,1)=zmat(i,idz)
zmat(i,idz)=temp
   end do
end if
!
!     Finally, update the matrix BMAT.
!
do j=1,n
   jp=npt+j
   w(jp)=bmat(knew,j)
   tempa=(alpha*vlag(jp)-tau*w(jp))/denom
   tempb=(-beta*w(jp)-tau*vlag(jp))/denom
   do i=1,jp
bmat(i,j)=bmat(i,j)+tempa*vlag(i)+tempb*w(i)
if (i > npt) bmat(jp,i-npt)=bmat(i,j)
   end do
end do
return
end
