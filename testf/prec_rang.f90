program prec_and_range
real(kind(0.0)) :: x
real(kind(0.0D0)) :: y(2)
integer(kind(0)) :: z

print *, precision(x), range(x)
print *, precision(y), range(y)
print *, range(z)
end program prec_and_range
