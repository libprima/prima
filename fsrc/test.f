      program test
          implicit none
          character(len = 1000) :: str
          character(len = 6) :: solver = 'newuoa' 
          character(len = 100) :: mssg 
          mssg='the objective function has been evaluated maxfun times.' 
          write(str, '(1X, 4A)')                                        &
     &     'Return from ', solver, ' because ', trim(mssg)
          print *, str          
          print *, trim(str)          
      end program test
