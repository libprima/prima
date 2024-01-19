
set(CMAKE_C_STANDARD 99)

include_directories(include)

if (WIN32 AND BUILD_SHARED_LIBS)
  add_compile_definitions(PRIMAC_EXPORTS)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    get_filename_component(Fortran_COMPILER_PATH ${CMAKE_Fortran_COMPILER} DIRECTORY)
    link_directories(${Fortran_COMPILER_PATH}/../../compiler/lib/intel64) # works for oneAPI 2023
  endif ()
  get_property(PRIMAF_MODULE_DIR TARGET primaf PROPERTY Fortran_MODULE_DIRECTORY)
  set(MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/../fortran/${PRIMAF_MODULE_DIR})
  if (CMAKE_CONFIGURATION_TYPES)
    set(MOD_DIR ${MOD_DIR}/${CMAKE_CFG_INTDIR})
  endif ()
  include_directories(${MOD_DIR})
endif ()

add_library(primac_f OBJECT cintrf.f90 cobyla_c.f90 lincoa_c.f90 bobyqa_c.f90 newuoa_c.f90 uobyqa_c.f90)
add_library(primac_c OBJECT prima.c)
add_library(primac $<TARGET_OBJECTS:primac_f> $<TARGET_OBJECTS:primac_c>)
add_dependencies(primac_f primaf)
add_dependencies(primac_c primaf)

if (WIN32)
  set_target_properties(primac PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()


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
  add_test (NAME example_${name}_c COMMAND example_${name}_c_exe)

  # Within CI, we'd like to run with gdb so that if there's a segfault the logs will have a stacktrace we can use to investigate.
  # Of course this can be run locally as well if you define CI in your environment.
  if(NOT APPLE AND UNIX AND DEFINED ENV{CI})  # Apple security policy will not allow running gdb in CI
    add_test (NAME example_${name}_c_with_gdb COMMAND gdb -batch --command=${CMAKE_BINARY_DIR}/cmdfile.gdb example_${name}_c_exe)
  elseif(WIN32 AND DEFINED ENV{CI})
    # For Windows we need to provide the full path to the executable since it is installed to a different directory
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
