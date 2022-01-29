PROGRAM testtype
 USE, INTRINSIC :: iso_fortran_env, ONLY: integer_kinds
 IMPLICIT NONE
 INTEGER, PARAMETER :: ik = kind(1)
 INTEGER, PARAMETER :: ikk = minval(integer_kinds, mask=integer_kinds/=ik)
 INTEGER(KIND=ik) i
 INTEGER(KIND=ikk) j
 i = 0
 j = 1
 k = 2
 WRITE (*, *) i, j, k
END PROGRAM testtype
