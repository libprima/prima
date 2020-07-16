      program test
          use iso_fortran_env
          implicit none
          print *, INT16, INT32, INT64
          print *, REAL16, REAL32, REAL64, REAL128
          print *, selected_real_kind(30, 291)
          print *, range(0.0_REAL128)
          print *, range(0.0_REAL64)
          print *, range(0.0_REAL32)
      end program test
