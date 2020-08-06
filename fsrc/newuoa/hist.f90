module hist

updatehist(ihist, )

! Inputs
integer(IK) :: ihist
integer(IK) :: nf
real(RP) :: f
real(RP) :: x(:)


! In-outputs
real(RP) :: fhist(:)
real(RP) :: xhist(:, :)

k = mod(nf - 1, size(fhist)) + 1
fhist(k) = f
xhist(:, k) = x

end if


end module hist
