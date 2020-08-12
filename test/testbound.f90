program testbound
implicit none
real :: a(3)
a = (/1., 2., 3./)
print *, '4:3', a(4:3)
print *, (/a(4:3), a(1:3)/)
print *, '1:0', a(1:0)
print *, (/a(4:3), a(1:3), a(1:0)/)
end program testbound
