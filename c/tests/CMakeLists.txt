
macro (prima_add_c_test_multi name)
  add_executable(${name}_c_exe EXCLUDE_FROM_ALL ${name}.c)
  if (PRIMA_ENABLE_TESTING)
    set_target_properties (${name}_c_exe PROPERTIES EXCLUDE_FROM_ALL FALSE)
  endif ()
  if (WIN32)
    set_target_properties(${name}_c_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
  endif()
  target_link_libraries(${name}_c_exe PRIVATE primac)
  add_dependencies(tests ${name}_c_exe)
  # Outside of CI we don't want to force people to run tests with gdb, so we test the executables by themselves.
  # We want these to run in CI as well, because sometimes running with gdb masks an error, so we set them up
  # before we set up the tests for CI
  add_test(NAME bobyqa_${name}_c COMMAND ${name}_c_exe bobyqa)
  add_test(NAME cobyla_${name}_c COMMAND ${name}_c_exe cobyla)
  add_test(NAME lincoa_${name}_c COMMAND ${name}_c_exe lincoa)
  add_test(NAME newuoa_${name}_c COMMAND ${name}_c_exe newuoa)
  add_test(NAME uobyqa_${name}_c COMMAND ${name}_c_exe uobyqa)

  # Within CI, we'd like to run with gdb so that if there's a segfault the logs will have a stacktrace we can use to investigate.
  # Of course this can be run locally as well if you define CI in your environment.
  if(NOT APPLE AND UNIX AND DEFINED ENV{CI})  # Apple security policy will not allow running gdb in CI
    add_test(NAME bobyqa_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${name}_c_exe bobyqa)
    add_test(NAME cobyla_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${name}_c_exe cobyla)
    add_test(NAME lincoa_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${name}_c_exe lincoa)
    add_test(NAME newuoa_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${name}_c_exe newuoa)
    add_test(NAME uobyqa_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${name}_c_exe uobyqa)
  elseif(WIN32 AND DEFINED ENV{CI})
    # For Windows we need to provide the full path to the executable since it is installed to a different directory
    add_test(NAME bobyqa_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${name}_c_exe.exe bobyqa)
    add_test(NAME cobyla_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${name}_c_exe.exe cobyla)
    add_test(NAME lincoa_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${name}_c_exe.exe lincoa)
    add_test(NAME newuoa_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${name}_c_exe.exe newuoa)
    add_test(NAME uobyqa_${name}_c_with_gdb COMMAND gdb  -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb --args ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/${name}_c_exe.exe uobyqa)
  endif()

endmacro ()

prima_add_c_test_multi(stress)
prima_add_c_test_multi(data)
