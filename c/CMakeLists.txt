
set(CMAKE_C_STANDARD 99)

add_library (primac cintrf.f90 cobyla_c.f90 lincoa_c.f90 bobyqa_c.f90 newuoa_c.f90 uobyqa_c.f90 prima.c)
if (WIN32)
  set_target_properties(primac PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()

# The newer intel compiler has issues when trying to capture a Fortran procedure pointer in a closure.
# We can fix this by appending norecursion to our compile options when we compile the C interface.
# More details may be found here: https://fortran-lang.discourse.group/t/strange-issue-with-ifx-compiler-and-assume-recursion/7013
if (CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  if (WIN32)
    # The $<$<COMPILE_LANGUAGE:Fortran> part is so that this only applies when compiling Fortran
    # files and *not* when compiling prima.c (because 'assume' is not a valid compiler flag when
    # compiling C).
    target_compile_options(primac PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:/assume:norecursion>)
  else ()
    target_compile_options(primac PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-assume norecursion>)
  endif ()
endif ()

target_include_directories (primac PUBLIC
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_link_libraries (primac PUBLIC primaf) # must be PUBLIC for precision macros
set_target_properties(primac PROPERTIES POSITION_INDEPENDENT_CODE ON C_STANDARD 99)

if (NOT BUILD_SHARED_LIBS)
  target_compile_definitions(primac PUBLIC PRIMAC_STATIC)
  target_link_libraries (primac INTERFACE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif ()

# Export symbols on Windows. See more comments in fortran/CMakeLists.txt.
if (WIN32 AND BUILD_SHARED_LIBS)
  target_sources(primac PRIVATE primac.def)
endif ()

install (FILES include/prima/prima.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prima)

macro (prima_add_c_test name)
  add_executable (example_${name}_c_exe EXCLUDE_FROM_ALL examples/${name}/${name}_example.c)
  target_link_libraries (example_${name}_c_exe PRIVATE primac)
  target_include_directories (example_${name}_c_exe PRIVATE ${CMAKE_SOURCE_DIR}/c/include)
  if (PRIMA_ENABLE_EXAMPLES)
    set_target_properties (example_${name}_c_exe PROPERTIES EXCLUDE_FROM_ALL FALSE)
  endif ()
  if (WIN32)
    set_target_properties(example_${name}_c_exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
  endif()

  # Outside of CI we don't want to force people to run examples with gdb, so we test the executables by themselves.
  # We want these to run in CI as well, because sometimes running with gdb masks an error, so we set them up
  # before we set up the examples for CI
  if(APPLE OR UNIX)
    add_test (NAME example_${name}_c COMMAND example_${name}_c_exe)
  elseif(WIN32)
    add_test (NAME example_${name}_c COMMAND ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/example_${name}_c_exe.exe)
  endif()

  # Within CI, we'd like to run with gdb so that if there's a segfault the logs will have a stacktrace we can use to investigate.
  # Of course this can be run locally as well if you define CI in your environment.
  if(NOT APPLE AND UNIX AND DEFINED ENV{CI})  # Apple security policy will not allow running gdb in CI
    add_test (NAME example_${name}_c_with_gdb COMMAND gdb -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb example_${name}_c_exe)
  elseif(WIN32 AND DEFINED ENV{CI})
    add_test (NAME example_${name}_c_with_gdb COMMAND gdb -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/example_${name}_c_exe.exe)
  endif()

  add_dependencies(examples example_${name}_c_exe)
endmacro ()

prima_add_c_test (cobyla)
prima_add_c_test (bobyqa)
prima_add_c_test (newuoa)
prima_add_c_test (uobyqa)
prima_add_c_test (lincoa)

add_subdirectory(tests)
