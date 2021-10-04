program testrand
implicit none
integer :: n
integer, allocatable :: seed0(:), seed1(:), seed2(:)
real :: r(3)

! Initialize the seed, and get its size.
call random_seed()
call random_seed(size=n)
allocate (seed0(n), seed1(n), seed2(n))

! Get the current seed.
call random_seed(get=seed0)

! Modify the seed, and put it as the new seed.
seed1 = max(seed0 - 1, 1)
call random_seed(put=seed1)

! Get the new seed.
call random_seed(get=seed2)

call random_seed(put=seed1)
call random_number(r)
write (*, *) r

call random_seed(put=seed2)
call random_number(r)
write (*, *) r

! Check whether the new seed is the one we put.
print *, maxval(abs(seed2 - seed1))

end program testrand
