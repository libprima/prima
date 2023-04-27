delsq = delta * delta
scaden = ZERO
biglsq = ZERO
knew = 0
do k = 1, npt
if (k == kopt) cycle
hdiag = ZERO
do jj = 1, nptm
hdiag = hdiag + zmat(k, jj)**2
end do
den = hdiag * beta + vlag(k)**2
distsq = ZERO
do j = 1, n
distsq = distsq + (xpt(j, k) - xopt(j))**2
end do
temp = max(ONE, (distsq / delsq)**2)
if (temp * den > scaden) then
scaden = temp * den
knew = k
denom = den
end if
if (temp * vlag(k)**2 > biglsq) biglsq = temp * vlag(k)**2
end do
if (.not. scaden > HALF * biglsq) then
if (nf > nresc) goto 190
!info = 4
info = DAMAGING_ROUNDING
goto 720
end if
