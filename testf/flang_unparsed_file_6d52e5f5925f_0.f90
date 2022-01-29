PROGRAM test
 USE :: consts_mod
 USE :: linalg_mod
 IMPLICIT NONE
 INTEGER, PARAMETER :: n = 5000
 REAL(KIND=rp) g(2,2)
 REAL(KIND=rp), PARAMETER :: a = sqrt(huge(0.0_rp))
 REAL(KIND=rp), PARAMETER :: x(2) = [a/pi, a/(pi*1.0e2_rp)]
 REAL(KIND=rp), PARAMETER :: y(2) = [x(2), x(1)]
 REAL(KIND=rp) r
 r = sqrt(sum(x**2))
 g = planerot(x)
 WRITE (*, *) "-->", maxval(abs(matprod(transpose(g), g)-eye(2))), norm(matpro&
 &d(g, x)-[r, zero])/max(r, eps)
 g = planerot(y)
 WRITE (*, *) "-->", maxval(abs(matprod(transpose(g), g)-eye(2))), norm(matpro&
 &d(g, y)-[r, zero])/max(r, eps)
END PROGRAM test
