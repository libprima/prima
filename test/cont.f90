program cont
    use iso_fortran_env, only: INT64, REAL64
    print *, real(0_INT64, REAL64)
    print *, real(1000_INT64, REAL64)
    print *, real(huge(0_INT64), REAL64)      
end program cont
