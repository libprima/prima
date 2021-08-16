program testmin
real :: a = 0
real :: b = 0
real :: c(2)
logical :: m(2) = [.true., .false.]
c = [a / a, b]
write (*, *) a, b
write (*, *) min(a, b)
write (*, *) minval([a, b])
write (*, *) min(b, b)
write (*, *) min(a, minval([b]))
write (*, *) min(a, minval(c, mask=m))
write (*, *) a, minval(c, mask=m)
end program testmin
