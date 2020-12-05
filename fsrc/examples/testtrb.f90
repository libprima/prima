program testtrb
!      use ifcore, only : tracebackqq
!    call tracebackqq(string='stop', user_exit_code=-1)
!  CALL TRACEBACKQQ(STRING="Bad value for TEMP",USER_EXIT_CODE=123)
    call backtrace()
    call abort
    print '(/1A/)', 'Finish.'
end program testtrb
