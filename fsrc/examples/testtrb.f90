program testtrb
      use ifcore, only : tracebackqq
    call tracebackqq(string='stop')
!  CALL TRACEBACKQQ(STRING="Bad value for TEMP",USER_EXIT_CODE=123)
!    call backtrace()
    print '(/1A/)', 'Finish.'
end program testtrb
